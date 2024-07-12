Skimming Workshop: Environment Setup
================
- <a href="#setting-up-the-pipeline"
  id="toc-setting-up-the-pipeline">Setting up the pipeline</a>

## Step 1: Detect Existing Conda/Miniforge/Miniconda Installation

Before proceeding with the installation of new software, check if there is any existing installation of Conda, Miniforge, or Miniconda on your system with the following command:

```bash
cat ~/.bash_profile
#export PATH="/Library/Python/3.9/lib/python/site-packages/radian"
cat ~/.condarc
#cat: /Users/sjr729/.condarc: No such file or directory
grep -i "conda" ~/.zshrc
#grep: /Users/sjr729/.zshrc: No such file or directory
cat ~/.zshrc
#cat: /Users/sjr729/.zshrc: No such file or directory
```

## Step 2: Install miniforge

1.  Install
    [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html)
    or conda
    (<https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html>)
    I'm doing a fresh install of mamba with
    miniforge (<https://github.com/conda-forge/miniforge>) because mamba is faster than conda. 
2. 
```bash
#wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
#bash Miniforge3-$(uname)-$(uname -m).sh
```
```html
#installation finished.
#Do you wish to update your shell profile to automatically initialize conda?
#This will activate conda on startup and change the command prompt when activated.
#If you'd prefer that conda's base environment not be activated on startup,
#run the following command when conda is activated:
#conda config --set auto_activate_base false
#You can undo this by running `conda init --reverse $SHELL`? [yes|no] yes
```
3. Check mamba/conda --version
```bash
#bash-3.2$ mamba --version
#mamba 1.5.8
#conda 24.3.0
```
4. do bash install.sh and name the conda env
```bash
cd ~/Desktop
git clone https://github.com/echarvel3/skimming_scripts-echarvel.git
bash ./install.sh
#Setting up Conda environment...
#Please enter:
#Name of New Conda Environment: -> workshop_test
```

```bash
==> Autoremoving 1 unneeded formula:
jsoncpp
Uninstalling /opt/homebrew/Cellar/jsoncpp/1.9.5... (18 files, 324.1KB)
Pruned 0 symbolic links and 4 directories from /opt/homebrew
==> Caveats
==> libomp
libomp is keg-only, which means it was not symlinked into /opt/homebrew,
because it can override GCC headers and result in broken builds.

For compilers to find libomp you may need to set:
  export LDFLAGS="-L/opt/homebrew/opt/libomp/lib"
  export CPPFLAGS="-I/opt/homebrew/opt/libomp/include"
==> llvm
To use the bundled libc++ please add the following LDFLAGS:
  LDFLAGS="-L/opt/homebrew/opt/llvm/lib/c++ -Wl,-rpath,/opt/homebrew/opt/llvm/lib/c++"

llvm is keg-only, which means it was not symlinked into /opt/homebrew,
because macOS already provides this software and installing another version in
parallel can cause all kinds of trouble.

If you need to have llvm first in your PATH, run:
  echo 'export PATH="/opt/homebrew/opt/llvm/bin:$PATH"' >> ~/.zshrc

For compilers to find llvm you may need to set:
  export LDFLAGS="-L/opt/homebrew/opt/llvm/lib"
  export CPPFLAGS="-I/opt/homebrew/opt/llvm/include"
g++-14-14 -Wno-unused-result -Wno-unused-command-line-argument -g -std=c++11 -O3 -fopenmp   -lm -lz -lstdc++ -lcurl -c src/io.cpp -o build/io.o
make: g++-14-14: No such file or directory
make: *** [build/io.o] Error 1
Downloading Bacterial Decontamination Library...
```
5. MANUAL INSTALLATION from scratch 
```bash
ENV_NAME="workshop_test"
#conda create --name $ENV_NAME --file ./Obsolete/respect-spec-file.txt
conda create --name $ENV_NAME python=3.7
conda activate $ENV_NAME
```
```
(base) bash-3.2$ conda create --name $ENV_NAME python=3.7
Channels:
 - conda-forge
Platform: osx-arm64
Collecting package metadata (repodata.json): done
Solving environment: failed

PackagesNotFoundError: The following packages are not available from current channels:

  - python=3.7*

Current channels:

  - https://conda.anaconda.org/conda-forge

To search for alternate channels that may provide the conda package you're
looking for, navigate to

    https://anaconda.org

and use the search bar at the top of the page.
```
