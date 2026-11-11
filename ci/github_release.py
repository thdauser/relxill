#!/usr/bin/env python3
"""Publish relxill on GitHub (https://github.com/thdauser/relxill), run by the GitLab CI (see .gitlab-ci.yml).

Two modes, both run only after build and tests passed:

  ci/github_release.py --push-master
      CI job "github-master", for every push to master: pushes the commit to GitHub master, if that is a
      fast-forward. GitHub master follows the development, and the documentation on GitHub Pages
      (.github/workflows/docs.yml) is rebuilt.

  ci/github_release.py vx.y.z
      CI job "github", for every tag vx.y.z:
        1. pushes the tag to GitHub
        2. creates (or updates) the GitHub release "version x.y.z", with the CHANGELOG.txt entry as text;
           it is marked as latest release unless a higher version was already released
        3. attaches the model tarball relxill_model_vx.y.z.tgz (replacing an existing one of the same name)

Both can be repeated safely (e.g., re-running a CI job) and run by hand from the top-level directory, e.g.,
for a tag made before the CI job existed:

    ci/github_release.py v2.10.0 --tarball path/to/relxill_model_v2.10.0.tgz --dry-run   # no token needed
    read -rs GITHUB_TOKEN && export GITHUB_TOKEN
    ci/github_release.py v2.10.0 --tarball path/to/relxill_model_v2.10.0.tgz

GITHUB_TOKEN: fine-grained token for the repository with the permissions "Contents" and "Workflows"
(read and write); the latter is needed to push changes of .github/workflows/.
Only the Python standard library is used.
"""

import argparse
import json
import os
import re
import subprocess
import sys
import tarfile
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

DEFAULT_REPO = "thdauser/relxill"
API = "https://api.github.com"
WEBPAGE = "https://www.sternwarte.uni-erlangen.de/research/relxill/"


def fail(msg):
    sys.exit(f"ERROR: {msg}")


def git(*args, check=True):
    return subprocess.run(["git", *args], check=check, capture_output=True, text=True)


def version_key(version):
    """(2, 10, 0) for "v2.10.0"; "v2.4" gives (2, 4, 0)."""
    parts = [int(p) for p in version.lstrip("v").split(".")]
    return tuple(parts + [0] * (3 - len(parts)))


# ---------------------------------------------------------------- release text

def changelog_entry(path, version):
    """Items of the entry "[x.y.z] - YYYY-MM-DD" in CHANGELOG.txt (see RELEASE.md), or None."""
    lines = Path(path).read_text().splitlines()
    head = re.compile(rf"^\[v?{re.escape(version)}\]")
    for start, line in enumerate(lines):
        if head.match(line):
            break
    else:
        return None
    items = []
    for line in lines[start + 1:]:
        if line.startswith("["):
            break
        text = line.strip()
        if text.startswith("- "):
            items.append(text[2:].strip())
        elif text and items:
            items[-1] += " " + text
    return items or None


def release_body(items, version):
    """Release text; None if CHANGELOG.txt has no entry (an existing release then keeps its text)."""
    if items is None:
        print(f"WARNING: no entry [{version}] with items in CHANGELOG.txt")
        return None
    if len(items) == 1:
        text = items[0]
    else:
        text = "\n".join(f"- {item}" for item in items)
    return (f"{text}\n\n"
            f"Download `relxill_model_v{version}.tgz` below and see the README for the installation. "
            f"The reflection tables are available from the [relxill webpage]({WEBPAGE}).\n")


# ---------------------------------------------------------------- GitHub API

class GitHub:
    def __init__(self, repo, token):
        self.repo, self.token = repo, token

    def request(self, method, url, data=None, content_type="application/json"):
        if not url.startswith("https://"):
            url = f"{API}/repos/{self.repo}{url}"
        headers = {"Accept": "application/vnd.github+json",
                   "X-GitHub-Api-Version": "2022-11-28",
                   "User-Agent": "relxill-release-script"}
        if self.token:
            headers["Authorization"] = f"Bearer {self.token}"
        if isinstance(data, dict):
            data = json.dumps(data).encode()
        if data is not None:
            headers["Content-Type"] = content_type
        req = urllib.request.Request(url, data=data, headers=headers, method=method)
        try:
            with urllib.request.urlopen(req, timeout=120) as resp:
                body = resp.read()
        except urllib.error.HTTPError as err:
            if err.code == 404 and method == "GET":
                return None
            fail(f"GitHub API {method} {url}: {err.code} {err.reason}\n{err.read().decode(errors='replace')}")
        return json.loads(body) if body else {}

    def highest_released_version(self, exclude_tag):
        releases = self.request("GET", "/releases?per_page=100") or []
        versions = [version_key(r["tag_name"]) for r in releases
                    if r["tag_name"] != exclude_tag and not r["draft"] and not r["prerelease"]
                    and re.fullmatch(r"v\d+(\.\d+){1,2}", r["tag_name"])]
        return max(versions, default=None)


# ---------------------------------------------------------------- git push

def git_push(repo, token, refspecs, dry_run):
    print(f"push to https://github.com/{repo}.git: {' '.join(refspecs)}")
    if dry_run:
        return
    auth_url = f"https://x-access-token:{token}@github.com/{repo}.git"
    result = git("push", auth_url, *refspecs, check=False)
    if result.returncode != 0:
        fail("git push to GitHub failed:\n" + (result.stdout + result.stderr).replace(token, "***"))


def push_master(repo, token, sha, dry_run):
    """Move GitHub master to `sha` if that is a fast-forward."""
    remote = git("ls-remote", f"https://github.com/{repo}.git", "refs/heads/master").stdout.split()
    gh_master = remote[0] if remote else None

    def is_ancestor(a, b):
        return git("merge-base", "--is-ancestor", a, b, check=False).returncode == 0

    if gh_master == sha:
        print(f"GitHub master is already at {sha[:9]}")
        return
    if gh_master is not None:
        if git("cat-file", "-e", f"{gh_master}^{{commit}}", check=False).returncode != 0:
            fail(f"GitHub master ({gh_master[:9]}) is not in this repository: something was committed directly "
                 f"on GitHub. Bring these changes into the GitLab repository first (see RELEASE.md).")
        if is_ancestor(sha, gh_master):
            print(f"GitHub master ({gh_master[:9]}) is already newer than {sha[:9]}, nothing to do")
            return
        if not is_ancestor(gh_master, sha):
            fail(f"GitHub master ({gh_master[:9]}) and {sha[:9]} have diverged (history rewritten?), "
                 f"not pushing. Check and fix GitHub master by hand.")
    git_push(repo, token, [f"{sha}:refs/heads/master"], dry_run)


# ---------------------------------------------------------------- release

def publish_release(gh, tag, version, body, tarball, dry_run):
    highest = gh.highest_released_version(exclude_tag=tag)
    latest = highest is None or version_key(version) >= highest
    fields = {"tag_name": tag, "name": f"version {version}",
              "draft": False, "prerelease": False, "make_latest": "true" if latest else "false"}
    release = gh.request("GET", f"/releases/tags/{urllib.parse.quote(tag)}")
    if body is None and not release:
        body = "See CHANGELOG.txt for the changes in this version.\n"
    if body is not None:
        fields["body"] = body
    action = "update" if release else "create"
    print(f"{action} GitHub release '{fields['name']}' (latest: {fields['make_latest']}), text:\n")
    print("    " + body.replace("\n", "\n    ") if body is not None else "    (unchanged)\n")
    if dry_run:
        print(f"upload asset {tarball.name} ({tarball.stat().st_size} bytes)")
        return None

    if release:
        release = gh.request("PATCH", f"/releases/{release['id']}", fields)
    else:
        release = gh.request("POST", "/releases", fields)

    for asset in release.get("assets", []):
        if asset["name"] == tarball.name:
            print(f"replace existing asset {asset['name']}")
            gh.request("DELETE", f"/releases/assets/{asset['id']}")
    upload_url = release["upload_url"].split("{")[0] + "?" + urllib.parse.urlencode({"name": tarball.name})
    gh.request("POST", upload_url, tarball.read_bytes(), content_type="application/gzip")
    print(f"upload asset {tarball.name} ({tarball.stat().st_size} bytes)")
    return release["html_url"]


def check_tarball(path, version):
    if not path.is_file():
        fail(f"tarball not found: {path}")
    with tarfile.open(path) as tar:
        try:
            config = tar.extractfile("config.h").read().decode()
        except KeyError:
            fail(f"{path} does not contain config.h, is it a relxill model tarball?")
    m = re.search(r'#define\s+PROJECT_VER\s+"([^"]+)"', config)
    if not m or m.group(1) != version:
        fail(f"{path} is version {m.group(1) if m else '?'}, not {version}")


# ---------------------------------------------------------------- main

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("tag", nargs="?", default=os.environ.get("CI_COMMIT_TAG"),
                        help="release tag vx.y.z (default: $CI_COMMIT_TAG)")
    parser.add_argument("--push-master", action="store_true",
                        help="push the current commit ($CI_COMMIT_SHA or HEAD) to GitHub master instead")
    parser.add_argument("--tarball", type=Path, help="model tarball (default: relxill_model_vx.y.z.tgz)")
    parser.add_argument("--changelog", type=Path, default=Path("CHANGELOG.txt"))
    parser.add_argument("--repo", default=DEFAULT_REPO, help="GitHub repository (default: %(default)s)")
    parser.add_argument("--dry-run", action="store_true", help="only show what would be done")
    args = parser.parse_args()

    token = os.environ.get("GITHUB_TOKEN", "")
    if not token and not args.dry_run:
        fail("GITHUB_TOKEN is not set.\n"
             "  by hand: run 'read -rs GITHUB_TOKEN && export GITHUB_TOKEN' in the same shell and paste the token\n"
             "  in GitLab CI: Settings > CI/CD > Variables (masked, protected); the branch master and the tags v*\n"
             "  must be protected for the variable to be available")

    if args.push_master:
        sha = os.environ.get("CI_COMMIT_SHA") or git("rev-parse", "HEAD").stdout.strip()
        push_master(args.repo, token, sha, args.dry_run)
        print("\n*** dry run, nothing was changed ***" if args.dry_run else "\ndone")
        return

    if not args.tag or not re.fullmatch(r"v\d+\.\d+\.\d+", args.tag):
        fail(f"need a tag of the form vx.y.z (or --push-master), got '{args.tag}'")
    tag, version = args.tag, args.tag[1:]

    if not git("rev-parse", "--verify", f"refs/tags/{tag}^{{commit}}", check=False).stdout.strip():
        fail(f"tag {tag} not found in this repository")

    tarball = args.tarball or Path(f"relxill_model_v{version}.tgz")
    check_tarball(tarball, version)
    body = release_body(changelog_entry(args.changelog, version), version)

    git_push(args.repo, token, [f"refs/tags/{tag}"], args.dry_run)
    url = publish_release(GitHub(args.repo, token), tag, version, body, tarball, args.dry_run)

    print("\n*** dry run, nothing was changed ***" if args.dry_run else f"\nreleased: {url}")


if __name__ == "__main__":
    main()
