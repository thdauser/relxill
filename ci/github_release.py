#!/usr/bin/env python3
"""Publish a relxill release on GitHub (https://github.com/thdauser/relxill).

Run by the GitLab CI job "github" for every tag vx.y.z (see .gitlab-ci.yml), after build and tests passed.
It can also be run by hand from the top-level directory, e.g., for a tag made before this job existed:

    GITHUB_TOKEN=... ci/github_release.py v2.10.0 --tarball path/to/relxill_model_v2.10.0.tgz
    ci/github_release.py v2.10.0 --tarball ... --dry-run      # no token needed, changes nothing

Steps (all of them can be repeated safely, e.g., when re-running the CI job):
  1. push the tag to GitHub, and move GitHub master to the tagged commit if that is a fast-forward
     (GitHub master always shows the latest release; for a tag on an older branch only the tag is pushed)
  2. create (or update) the GitHub release "version x.y.z", with the CHANGELOG.txt entry as text
  3. attach the model tarball relxill_model_vx.y.z.tgz (replacing an existing one of the same name)

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


def push_to_github(repo, token, tag, sha, dry_run):
    """Push the tag; move GitHub master to the tag if it is a fast-forward. Return True if master moved."""
    public_url = f"https://github.com/{repo}.git"
    remote = git("ls-remote", public_url, "refs/heads/master").stdout.split()
    gh_master = remote[0] if remote else None

    if gh_master == sha:
        update_master = False
        print(f"GitHub master is already at {sha[:9]}")
    elif gh_master is None:
        update_master = True
    elif git("cat-file", "-e", f"{gh_master}^{{commit}}", check=False).returncode != 0:
        update_master = False
        print(f"WARNING: GitHub master ({gh_master[:9]}) is not in this repository (something was committed "
              f"directly on GitHub?). Only pushing the tag; bring the GitHub changes into GitLab first.")
    else:
        update_master = git("merge-base", "--is-ancestor", gh_master, sha, check=False).returncode == 0
        if not update_master:
            print(f"{tag} is not newer than GitHub master ({gh_master[:9]}), only pushing the tag")

    refspecs = [f"refs/tags/{tag}"] + ([f"{sha}:refs/heads/master"] if update_master else [])
    print(f"push to {public_url}: {' '.join(refspecs)}")
    if not dry_run:
        auth_url = f"https://x-access-token:{token}@github.com/{repo}.git"
        result = git("push", auth_url, *refspecs, check=False)
        if result.returncode != 0:
            fail("git push to GitHub failed:\n" + (result.stdout + result.stderr).replace(token, "***"))
    return update_master or gh_master == sha


def publish_release(gh, tag, version, body, latest, tarball, dry_run):
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


# ---------------------------------------------------------------- main

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


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("tag", nargs="?", default=os.environ.get("CI_COMMIT_TAG"),
                        help="tag vx.y.z (default: $CI_COMMIT_TAG)")
    parser.add_argument("--tarball", type=Path, help="model tarball (default: relxill_model_vx.y.z.tgz)")
    parser.add_argument("--changelog", type=Path, default=Path("CHANGELOG.txt"))
    parser.add_argument("--repo", default=DEFAULT_REPO, help="GitHub repository (default: %(default)s)")
    parser.add_argument("--dry-run", action="store_true", help="only show what would be done")
    args = parser.parse_args()

    if not args.tag or not re.fullmatch(r"v\d+\.\d+\.\d+", args.tag):
        fail(f"need a tag of the form vx.y.z, got '{args.tag}'")
    tag, version = args.tag, args.tag[1:]

    token = os.environ.get("GITHUB_TOKEN", "")
    if not token and not args.dry_run:
        fail("GITHUB_TOKEN is not set (in GitLab: Settings > CI/CD > Variables, masked and protected, "
             "and the tag v* must be protected for the variable to be available)")

    sha = git("rev-parse", "--verify", f"refs/tags/{tag}^{{commit}}", check=False).stdout.strip()
    if not sha:
        fail(f"tag {tag} not found in this repository")

    tarball = args.tarball or Path(f"relxill_model_v{version}.tgz")
    check_tarball(tarball, version)
    body = release_body(changelog_entry(args.changelog, version), version)

    latest = push_to_github(args.repo, token, tag, sha, args.dry_run)
    url = publish_release(GitHub(args.repo, token), tag, version, body, latest, tarball, args.dry_run)

    print("\n*** dry run, nothing was changed ***" if args.dry_run else f"\nreleased: {url}")


if __name__ == "__main__":
    main()
