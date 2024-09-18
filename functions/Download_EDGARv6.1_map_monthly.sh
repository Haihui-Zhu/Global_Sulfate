#!/bin/bash

# This script downloads EDGARv6.1 air pollutant monthly map data

exit_with_error() {
    echo
    echo "Unable to Retrieve Data"
    echo
    echo $1
    echo
    echo "Please check the URL or network settings."
    echo
    exit 1
}

fetch_urls() {
  if command -v wget >/dev/null 2>&1; then
      echo
      echo "Using wget to download files."
      echo
      while IFS= read -r line; do
        # Get everything after the last '/'
        filename="${line##*/}"

        # Strip everything after '?'
        stripped_query_params="${filename%%\?*}"

        # Determine the directory based on the URL
        dir="/nobackup/hzhu3/ExtData/EDGARv6.1"

        # Full path to where the file would be stored
        full_path="$dir/$stripped_query_params"

        # Check if the file already exists
        if [[ -f "$full_path" ]]; then
            echo "File already exists, skipping: $full_path"
            continue
        fi

        # Base name of the file (without extension)
        base_name="${stripped_query_params%.*}"

        # Check if a directory with the base name exists
        if [[ -d "$dir/$base_name" ]]; then
            echo "Directory already exists, skipping: $dir/$base_name"
            continue
        fi

        # Download the file to the appropriate directory
        wget --output-document "$full_path" -- "$line" && echo || exit_with_error "Command failed with error. Please retrieve the data manually."
      done < "EDGAR_download_urls.txt"
  else
      exit_with_error "Error: wget not found. Please install wget."
  fi
}


fetch_urls
