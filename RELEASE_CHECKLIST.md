# PyPI Release Checklist

## How GitHub Hooks Up with PyPI

Your project uses **GitHub Actions** to automatically publish to PyPI when you create a GitHub release. Here's how it works:

1. **Workflow File**: `.github/workflows/build.yml` contains the publishing workflow
2. **Trigger**: The workflow runs automatically when you create a GitHub release
3. **Authentication**: Uses a PyPI API token stored as a GitHub secret (`PYPI_API_TOKEN`)
4. **Action**: Uses `pypa/gh-action-pypi-publish@release/v1` which is the official PyPI publishing action

## Current Configuration

- **Workflow**: `.github/workflows/build.yml`
- **Trigger**: GitHub releases (when you create a release on GitHub)
- **Secret Required**: `PYPI_API_TOKEN`

## Step 1: Verify PyPI API Token Secret is Set

1. Go to: https://github.com/hstern2/amsr/settings/secrets/actions
2. Check if `PYPI_API_TOKEN` exists in the list
3. If it doesn't exist, you'll need to create it (see Step 2)

## Step 2: Create PyPI API Token (if needed)

If `PYPI_API_TOKEN` doesn't exist, you need to create it:

1. **Log in to PyPI**: https://pypi.org/account/login/
2. **Go to Account Settings**: https://pypi.org/manage/account/
3. **Scroll to "API tokens"** section
4. **Click "Add API token"**
5. **Name it** (e.g., "GitHub Actions - amsr")
6. **Scope**: Select "Entire account" or "Project: amsr" (if you have project-specific tokens)
7. **Copy the token** (you'll only see it once!)
8. **Add to GitHub Secrets**:
   - Go to: https://github.com/hstern2/amsr/settings/secrets/actions
   - Click "New repository secret"
   - Name: `PYPI_API_TOKEN`
   - Value: Paste your PyPI API token
   - Click "Add secret"

## Step 3: Test the Workflow (Dry Run)

1. **Create a TestPyPI API token**:
   - Go to: https://test.pypi.org/manage/account/
   - Create an API token (same process as above)

2. **Add TestPyPI secret to GitHub**:
   - Name: `TEST_PYPI_API_TOKEN`
   - Value: Your TestPyPI API token

3. **Temporarily modify `.github/workflows/build.yml`**:
   - Uncomment line 37: `repository-url: https://test.pypi.org/legacy/`
   - Change line 35 to use `TEST_PYPI_API_TOKEN` instead of `PYPI_API_TOKEN`

4. **Create a test release**:
   - Create a new tag: `git tag v0.1.6-test`
   - Push the tag: `git push origin v0.1.6-test`
   - Create a GitHub release for that tag
   - Watch the workflow run at: https://github.com/hstern2/amsr/actions

5. **Verify on TestPyPI**: Check https://test.pypi.org/project/amsr/

6. **Revert the workflow changes** before real release

## Step 4: Prepare for Real Release

1. **update version**

   1. `api.rst`
   2. `pyproject.toml`
   3. `amsr/version.py`
   4. `docs/source/changelog.rst`
   5. `docs/source/conf.py`

2. **commit and push**

3. **create a GitHub release**:
   - Go to: https://github.com/hstern2/amsr/releases/new
   - Create a new tag (e.g., `v0.1.7`)
   - Add release notes
   - Click "Publish release"

4. **monitor the workflow**:
   - Watch: https://github.com/hstern2/amsr/actions
   - The workflow should automatically start and publish to PyPI

5. **verify on PyPI**:
   - Check: https://pypi.org/project/amsr/
   - Your new version should appear within a few minutes
