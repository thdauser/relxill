# Releasing a new relxill version

Checklist for maintainers. A release touches three places:

* `relxill-code` (this repository, on the internal GitLab): version number, change log, tag, model tarball
* https://github.com/thdauser/relxill: the public copy with the releases and the documentation,
  updated automatically by the GitLab CI (GitHub `master` follows GitLab `master`, see section 2)
* `relxill-webpage` (`git@serpens.sternwarte.uni-erlangen.de:relxill/relxill-webpage.git`):
  the webpage https://www.sternwarte.uni-erlangen.de/research/relxill/ with the downloads

The version scheme (x.y.z) is explained at the top of `CHANGELOG.txt`.


## 1. Change log and version number (relxill-code)

1. In `CHANGELOG.txt`, rename the `[next version]` section to the new version and the release date:

   ```
   [2.10.0] - 2026-10-20
    - first change, can continue
      on the next line
    - second change
   ```

   Keep exactly this header format (full x.y.z version, date on the header line, one ` - ` per item):
   the webpage script and the GitHub release (`ci/github_release.py`) read the entry from here and turn each
   item into one bullet point.
   Start a new empty `[next version]` section above it for the changes after the release.
2. Set the version in `project(relxill VERSION x.y.z)` in the top-level `CMakeLists.txt`.
3. Commit with the message `[x.y.z] version update`.


## 2. Tag: tarball and GitHub release (relxill-code)

```bash
git tag -a vx.y.z -m "version x.y.z"
git push origin master vx.y.z
```

On a tag, the GitLab CI (`.gitlab-ci.yml`) builds and tests the model. If the tests pass, two jobs run:

* `provide` copies `relxill_model_vx.y.z.tgz` (and the `dev` tarball) to `/home/dauser/work/relxill/relxillArchive/`,
  where the webpage script takes it from.
* `github` (only for tags of the form `vx.y.z`) runs `ci/github_release.py vx.y.z`, which pushes the tag to
  https://github.com/thdauser/relxill and creates the GitHub release "version x.y.z" with the `CHANGELOG.txt` entry
  as text and the tarball attached. It is marked as the latest release unless a higher version exists.

Check that both jobs succeeded. Always release the tarball built by the CI from the tag, not one built locally
from a later commit. Both jobs can be re-run from the GitLab pipeline page if something went wrong.

To publish a release by hand (e.g., for a tag made before the `github` job existed):

```bash
ci/github_release.py vx.y.z --tarball /home/dauser/work/relxill/relxillArchive/relxill_model_vx.y.z.tgz --dry-run
read -rs GITHUB_TOKEN && export GITHUB_TOKEN     # paste the token, it is not shown
ci/github_release.py vx.y.z --tarball /home/dauser/work/relxill/relxillArchive/relxill_model_vx.y.z.tgz
```

### GitHub master and the documentation

Independent of releases, the CI job `github-master` pushes every commit on GitLab `master` to GitHub `master`
(`ci/github_release.py --push-master`) once its tests passed. GitHub `master` therefore follows the development,
and each push rebuilds the Doxygen documentation on GitHub Pages (`.github/workflows/docs.yml`). Users who want a
stable version take a release (tag or tarball).

Only fast-forwards are pushed. Pull requests on GitHub are therefore not merged there: bring them into the GitLab
repository (e.g., add the fork as remote and merge or cherry-pick the branch). If GitHub `master` contains a commit
unknown to GitLab, the `github-master` job fails until that is resolved.


## 3. Update the webpage (relxill-webpage)

The script `release_version.py` in the webpage repository does the routine part. It expects
`relxill-code` and `relxill-webpage` to be checked out next to each other (otherwise use `--changelog`).

```bash
cd ../relxill-webpage
git checkout -b update_vx.y.z
./release_version.py x.y.z            # dry run: shows the diff of index.html
./release_version.py x.y.z --apply
```

It copies the tarball into `relxillArchive/` and links it, and in `index.html` it

* sets the version number in the banner,
* makes the new tarball the latest download and moves the previous one to "Older Model Versions",
* adds the change log entry from `CHANGELOG.txt` and keeps only the newest 4 entries visible.

Then edit by hand:

* the banner text after "New version x.y",
* the "What's new" list (only for releases worth highlighting),
* the wording of the new change log entry, if needed (e.g., make the main statement bold).

Run `./check_webpage.py` (local checks), check the page in a browser, commit, push, and merge. To deploy, run as user
`dauser` on the remeis cluster:

```bash
./update_webpage.sh          # dry run
./update_webpage.sh --copy   # copy to the web server directory, then check the live page
```

After copying, the live page is checked automatically (`./check_webpage.py --live`, can also be run by hand from
anywhere): it must be identical to the local `index.html`, and all its downloads must work.

More details are in the `README.md` of the webpage repository.


## One-time setup

For the CI jobs `github-master` and `github`:

1. On GitHub, create a fine-grained personal access token (Settings > Developer settings > Personal access tokens >
   Fine-grained tokens) for the repository `thdauser/relxill` only, with the permissions "Contents: Read and write"
   and "Workflows: Read and write" (GitHub rejects pushes that change `.github/workflows/` without the latter).
   Note its expiry date: the job fails once it expired.
2. In GitLab (relxill-code > Settings > CI/CD > Variables), add it as `GITHUB_TOKEN`, with "Masked" and "Protected".
3. In GitLab (Settings > Repository > Protected tags), protect the tags `v*`. Protected variables are only passed to
   pipelines of protected branches and tags. Also check that `master` is listed under Protected branches (the
   default branch usually is).
4. On GitHub (thdauser/relxill > Settings > Pages), set the source to "GitHub Actions" for the documentation.

The CI runner needs `python3` and HTTPS access to github.com.
