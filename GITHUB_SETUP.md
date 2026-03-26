# How to Push to GitHub

Follow these steps to upload your files to GitHub.

## Step 1: Create Private Repository on GitHub

1. Go to https://github.com/new
2. Fill in:
   - **Repository name:** `Graphcomm-Lite-Cancer-Analysis`
   - **Description:** "AI-driven cancer cell communication analysis pipeline"
   - **Privacy:** ✅ **Private**
   - **Don't** initialize README / .gitignore / license (we already have them)
3. Click **"Create repository"**
4. You'll see setup instructions. Copy your repository URL (looks like: `https://github.com/YOUR_USERNAME/Graphcomm-Lite-Cancer-Analysis.git`)

## Step 2: Initialize Git Locally & Push

Open PowerShell in your project folder and run these commands:

```powershell
# Navigate to your project
cd "C:\Users\hamna pc\OneDrive\Desktop\GraphComm_Lite_Project"

# Initialize git repository
git init

# Add all files
git add .

# Create initial commit
git commit -m "Initial commit: GraphComm-Lite pipeline with 9 modular components"

# Add remote (replace YOUR_USERNAME)
git remote add origin https://github.com/YOUR_USERNAME/Graphcomm-Lite-Cancer-Analysis.git

# Push to GitHub (main branch)
git branch -M main
git push -u origin main
```

## Step 3: Verify on GitHub

1. Go to https://github.com/YOUR_USERNAME/Graphcomm-Lite-Cancer-Analysis
2. Confirm all files are there:
   - ✅ main.py
   - ✅ data_preprocess.py
   - ✅ snrna_features.py
   - ✅ pathway_analysis.py
   - ✅ graph_model.py
   - ✅ train_predict.py
   - ✅ drug_integration.py
   - ✅ relapse_analysis.py
   - ✅ visualize.py
   - ✅ requirements.txt
   - ✅ README.md
   - ✅ .gitignore
   - ✅ SNRNA_INTEGRATION.md

## Step 4: Share Link with Professor

Your repository link is:
```
https://github.com/YOUR_USERNAME/Graphcomm-Lite-Cancer-Analysis
```

Send this to your professor. They can:
- ✅ View the code (it's private, only they can see it with the link if you add them)
- ✅ Clone it: `git clone https://github.com/YOUR_USERNAME/Graphcomm-Lite-Cancer-Analysis.git`

---

## For First-Time GitHub Users

### Option A: Using GitHub Desktop (Easiest)

1. Download: https://desktop.github.com/
2. Sign in with your GitHub account
3. Click "(+)" → "Create" 
4. Choose your project folder
5. Click "Publish repository" → Make it Private
6. Done!

### Option B: Generate Personal Access Token (For Command Line)

If `git push` asks for password:

1. Go to GitHub → Settings → Personal access tokens → Tokens (classic)
2. Click "Generate new token (classic)"
3. Give it `repo` permissions
4. Copy the token
5. When prompted for password in PowerShell, paste this token

---

## Common Issues

### "fatal: remote origin already exists"
```powershell
git remote remove origin
git remote add origin https://github.com/YOUR_USERNAME/Graphcomm-Lite-Cancer-Analysis.git
git push -u origin main
```

### "error: The following untracked working tree files would be overwritten"
```powershell
git clean -fd
git reset --hard
```

### "fatal: not a git repository"
Make sure you ran `git init` in the correct folder first.

---

**Need help?** Paste the exact error message and I can troubleshoot!
