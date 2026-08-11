#!/bin/bash

set -e

# Directory where all operations will be performed
WORKING_DIR="GENOMICON-Seq"

# Arguments
MODE="${1:-}"
VERSION="${2:-v1.2.2}"

# Accept version with or without "v"
# 1.2.2  -> v1.2.2
# v1.2.2 -> v1.2.2
VERSION="v${VERSION#v}"

# Accept mode with or without "--"
case "$MODE" in
    wes|--wes)
        MODE="--wes"
        ;;
    ampliseq|--ampliseq)
        MODE="--ampliseq"
        ;;
    "")
        MODE=""
        ;;
    *)
        echo "Unknown mode: $MODE"
        echo "Use: wes, --wes, ampliseq, or --ampliseq"
        exit 1
        ;;
esac

# Create and enter the working directory
mkdir -p "$WORKING_DIR"
cd "$WORKING_DIR"

# Docker image and Git repository details
DOCKER_IMAGE="mimsto86/genomicon-seq:${VERSION}"
GIT_REPO_URL="https://github.com/Rounge-lab/GENOMICON-Seq"
RAW_BASE_URL="https://raw.githubusercontent.com/Rounge-lab/GENOMICON-Seq/${VERSION}"

FOLDER_TO_CLONE=(
    "input_data_ampliseq"
    "input_data_wes"
)

FILES_TO_DOWNLOAD=(
    "${RAW_BASE_URL}/Snakefile_ampliseq"
    "${RAW_BASE_URL}/Snakefile_wes"
    "${RAW_BASE_URL}/config_ampliseq.yml"
    "${RAW_BASE_URL}/config_wes.yml"
)

echo "Using GENOMICON-Seq version: ${VERSION}"

# Pull matching Docker image
echo "Pulling Docker image: ${DOCKER_IMAGE}"
docker pull "${DOCKER_IMAGE}"

# Function to download a specific folder from the selected GitHub release
clone_and_move() {
    local folder=$1

    echo "Downloading ${folder} from GENOMICON-Seq ${VERSION}"

    rm -rf myrepo

    git init myrepo
    cd myrepo

    git remote add origin "${GIT_REPO_URL}"
    git config core.sparseCheckout true

    echo "${folder}/" > .git/info/sparse-checkout

    # Fetch selected release/tag
    git fetch --depth 1 origin \
        "refs/tags/${VERSION}:refs/tags/${VERSION}"

    git checkout "${VERSION}"

    cd ..

    mv "myrepo/${folder}" .
    rm -rf myrepo
}

# Function to download individual files
download_files() {
    local files=("$@")

    echo "Downloading files from GENOMICON-Seq ${VERSION}"

    for file in "${files[@]}"; do
        curl -fLO "$file"
    done
}

# Setup based on selected mode
case "$MODE" in

    --ampliseq)

        clone_and_move "input_data_ampliseq"

        download_files \
            "${RAW_BASE_URL}/Snakefile_ampliseq" \
            "${RAW_BASE_URL}/config_ampliseq.yml"

        mkdir -p SQL_database
        cd SQL_database

        echo "Downloading SQLite database file for ampliseq"

        wget \
            https://zenodo.org/records/12683302/files/HPV16REF.sqlite \
            -O HPV16REF.sqlite

        cd ..
        ;;

    --wes)

        clone_and_move "input_data_wes"

        download_files \
            "${RAW_BASE_URL}/Snakefile_wes" \
            "${RAW_BASE_URL}/config_wes.yml"

        mkdir -p SQL_database
        cd SQL_database

        echo "Downloading SQLite database file for WES"

        wget \
            https://zenodo.org/records/12683302/files/chr1.sqlite \
            -O chr1.sqlite

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

        wget \
            https://zenodo.org/records/12683302/files/chr1.sqlite \
            -O chr1.sqlite

        wget \
            https://zenodo.org/records/12683302/files/HPV16REF.sqlite \
            -O HPV16REF.sqlite

        cd ..
        ;;
esac

echo
echo "Setup complete."
echo "GENOMICON-Seq version: ${VERSION}"
echo "Docker image: ${DOCKER_IMAGE}"
