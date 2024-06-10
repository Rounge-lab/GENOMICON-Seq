#!/bin/bash

# Docker image and Git repository details
DOCKER_IMAGE="mimsto86/genomicon-seq:v1.0"
GIT_REPO_URL="https://github.com/Rounge-lab/GENOMICON-Seq"
FOLDER_TO_CLONE=("input_data_ampliseq" "input_data_wes")
FILES_TO_DOWNLOAD=("https://raw.githubusercontent.com/Rounge-lab/GENOMICON-Seq/main/Snakefile_ampliseq" "https://raw.githubusercontent.com/Rounge-lab/GENOMICON-Seq/main/Snakefile_wes" "https://raw.githubusercontent.com/Rounge-lab/GENOMICON-Seq/main/config_ampliseq.yml" "https://raw.githubusercontent.com/Rounge-lab/GENOMICON-Seq/main/config_wes.yml")

# Pull the Docker image
echo "Pulling Docker image: $DOCKER_IMAGE"
docker pull $DOCKER_IMAGE

# Clone a specific folder from the Git repository
echo "Cloning specific folder from Git repository: $GIT_REPO_URL"
git init myrepo
cd myrepo
git remote add -f origin $GIT_REPO_URL
git config core.sparseCheckout true
for folder in "${FOLDER_TO_CLONE[@]}"; do
    echo "${folder}/" >> .git/info/sparse-checkout
done
git checkout main

# Move the folder to the desired location and cleanup
cd ..
for folder in "${FOLDER_TO_CLONE[@]}"; do
    mv "myrepo/${folder}" .
done
rm -rf myrepo

# Download individual files
echo "Downloading individual files from Git repository"
for file in "${FILES_TO_DOWNLOAD[@]}"; do
    curl -O "$file"
    # or use wget if you prefer: wget "$file"
done

echo "Setup complete."
