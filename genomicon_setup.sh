#!/bin/bash

# Directory where all operations will be performed
WORKING_DIR="GENOMICON-Seq"

# Create and enter the working directory
mkdir -p "$WORKING_DIR"
cd "$WORKING_DIR"

# Docker image and Git repository details
DOCKER_IMAGE="mimsto86/genomicon-seq:v1.0"
GIT_REPO_URL="https://github.com/Rounge-lab/GENOMICON-Seq"
FOLDER_TO_CLONE=("input_data_ampliseq" "input_data_wes")
FILES_TO_DOWNLOAD=("https://raw.githubusercontent.com/Rounge-lab/GENOMICON-Seq/main/Snakefile_ampliseq" "https://raw.githubusercontent.com/Rounge-lab/GENOMICON-Seq/main/Snakefile_wes" "https://raw.githubusercontent.com/Rounge-lab/GENOMICON-Seq/main/config_ampliseq.yml" "https://raw.githubusercontent.com/Rounge-lab/GENOMICON-Seq/main/config_wes.yml")

# Argument parsing
MODE=$1

# Pull the Docker image
echo "Pulling Docker image: $DOCKER_IMAGE"
docker pull $DOCKER_IMAGE

# Function to clone and move specific folders
clone_and_move() {
    local folder=$1
    echo "Cloning specific folder from Git repository: $folder"
    git init myrepo
    cd myrepo
    git remote add -f origin $GIT_REPO_URL
    git config core.sparseCheckout true
    echo "${folder}/" >> .git/info/sparse-checkout
    git checkout main
    cd ..
    mv "myrepo/${folder}" .
    rm -rf myrepo
}

# Function to download specific files
download_files() {
    local files=("$@")
    echo "Downloading individual files from Git repository"
    for file in "${files[@]}"; do
        curl -O "$file"
    done
}

# Download and setup based on the mode
case $MODE in
    --ampliseq)
        clone_and_move "input_data_ampliseq"
        download_files "https://raw.githubusercontent.com/Rounge-lab/GENOMICON-Seq/main/Snakefile_ampliseq" "https://raw.githubusercontent.com/Rounge-lab/GENOMICON-Seq/main/config_ampliseq.yml"
        mkdir -p SQL_database
        cd SQL_database
        echo "Downloading SQLite database file for ampliseq"
        wget https://zenodo.org/records/12683302/files/HPV16REF.sqlite -O HPV16REF.sqlite
        cd ..
        ;;
    --wes)
        clone_and_move "input_data_wes"
        download_files "https://raw.githubusercontent.com/Rounge-lab/GENOMICON-Seq/main/Snakefile_wes" "https://raw.githubusercontent.com/Rounge-lab/GENOMICON-Seq/main/config_wes.yml"
        mkdir -p SQL_database
        cd SQL_database
        echo "Downloading SQLite database file for wes"
        wget https://zenodo.org/records/12683302/files/chr1.sqlite -O chr1.sqlite
        cd ..
        ;;
    *)
        # Default behavior: download everything
        for folder in "${FOLDER_TO_CLONE[@]}"; do
            clone_and_move "$folder"
        done
        download_files "${FILES_TO_DOWNLOAD[@]}"
        mkdir -p SQL_database
        cd SQL_database
        echo "Downloading SQLite database files"
        wget https://zenodo.org/records/12683302/files/chr1.sqlite -O chr1.sqlite
        wget https://zenodo.org/records/12683302/files/HPV16REF.sqlite -O HPV16REF.sqlite
        cd ..
        ;;
esac

echo "Setup complete."
