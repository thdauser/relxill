# Releasing a new relxill version

Checklist for maintainers. A release touches three places:

* `relxill-code` (this repository, on the internal GitLab): version number, change log, tag, model tarball
* https://github.com/thdauser/relxill: the public copy with the releases and the documentation,
  updated automatically by the GitLab CI (GitHub `master` follows GitLab `master`, see section 2)
* the webpage https://www.sternwarte.uni-erlangen.de/research/relxill/ with the downloads, kept in the internal
  repository `relxill-webpage`

The version scheme (x.y.z) is explained at the top of `CHANGELOG.txt`. The internal parts of a release (the release
archive, deploying the webpage, publishing a release by hand, and the one-time setup of the CI) are described in the
`README.md` of `relxill-webpage`.


## 1. Change log and version number (relxill-code)

`master` is protected on GitLab: direct pushes are rejected, so the version update goes in through a merge request
like any other change.

1. Start a branch from an up-to-date `master`:

   ```bash
   git switch master && git pull
   git switch -c release-x.y.z
   ```

2. In `CHANGELOG.txt`, rename the `[next version]` section to the new version and the release date:

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
3. Set the version in `project(relxill VERSION x.y.z)` in the top-level `CMakeLists.txt`.
4. Commit with the message `[x.y.z] version update`, push the branch, and merge it into `master` with a merge
   request:

   ```bash
   git push -u origin release-x.y.z
   ```


## 2. Tag: tarball and GitHub release (relxill-code)

Tag only after the merge request is merged, and tag the commit on `master`. GitLab may create new commits when
merging (a merge commit, or a rebase), so a tag set on the release branch before can end up on a commit that is not
on `master`. For the same reason, never push the tag together with a branch: if the branch is rejected, the tag is
still pushed and starts the release (a push is only all-or-nothing with `--atomic`).

```bash
git switch master && git pull
git log -1 --oneline                  # must be "[x.y.z] version update" or the merge of it
git tag -a vx.y.z -m "version x.y.z"
git push origin vx.y.z
```

On a tag, the GitLab CI (`.gitlab-ci.yml`) builds and tests the model. If the tests pass, two jobs run:

* `provide` copies `relxill_model_vx.y.z.tgz` (and the `dev` tarball) to the release archive, from where the webpage
  takes it (section 3).
* `github` (only for tags of the form `vx.y.z`) runs `ci/github_release.py vx.y.z`, which pushes the tag to
  https://github.com/thdauser/relxill and creates the GitHub release "version x.y.z" with the `CHANGELOG.txt` entry
  as text and the tarball attached. It is marked as the latest release unless a higher version exists.

Check that both jobs succeeded. Always release the tarball built by the CI from the tag, not one built locally
from a later commit. Both jobs can be re-run from the GitLab pipeline page if something went wrong (publishing a
release by hand is described in the `README.md` of `relxill-webpage`).

### GitHub master and the documentation

Independent of releases, the CI job `github-master` pushes every commit on GitLab `master` to GitHub `master`
(`ci/github_release.py --push-master`) once its tests passed. GitHub `master` therefore follows the development,
and each push rebuilds the Doxygen documentation on GitHub Pages (`.github/workflows/docs.yml`). Users who want a
stable version take a release (tag or tarball).

Only fast-forwards are pushed. Pull requests on GitHub are therefore not merged there: bring them into the GitLab
repository (e.g., add the fork as remote and merge or cherry-pick the branch). If GitHub `master` contains a commit
unknown to GitLab, the `github-master` job fails until that is resolved.


## 3. Update the webpage (relxill-webpage)

The webpage is updated with the scripts in `relxill-webpage`, which expect `relxill-code` to be checked out next to
it: `release_version.py` adds the new tarball and the change log entry to `index.html`, `check_webpage.py` checks the
page, and `update_webpage.sh` deploys it. The steps are described in the `README.md` of `relxill-webpage`.
