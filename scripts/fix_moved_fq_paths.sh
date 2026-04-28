#!/usr/bin/env bash

# 2026-02-06 JCF/ChatGPT..

# Purpose:
# Sample table contains paths to fq files, but over time things may get rearranged.
# Typically, I'll know the new parent directory of the files, so I just want to
# take the basename from columns 4 and 5 of the sample table, find them in a new
# directory, and write a new file with the corrected names.


# Exit immediately if:
#  - any command fails (-e)
#  - any variable is undefined (-u)
#  - any command in a pipeline fails (pipefail)
set -euo pipefail

# Check for correct number of command-line arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input.tsv> <new_root_dir> <output.tsv>"
    exit 1
fi

# Assign command-line arguments to variables
INPUT="$1"       # Original TSV file
NEWROOT="$2"     # New directory tree where files now live
OUTPUT="$3"      # Output TSV file with corrected paths

# Run awk to process the file
awk -F'\t' -v OFS='\t' -v root="$NEWROOT" '

# Function to extract the basename (filename only) from a path
function basename(path,   n,a) {
    n = split(path, a, "/")
    return a[n]
}

# Function to find exactly one matching file in the new root directory
# - Returns the full path if exactly one match is found
# - Returns "NOT_FOUND" if no match is found
# - Exits with error if more than one match is found
function find_unique(basename,   cmd, line, count, result) {
    cmd = "find \"" root "\" -type f -name \"" basename "\" 2>/dev/null"

    count = 0
    result = ""

    while ((cmd | getline line) > 0) {
        count++
        result = line

        # If more than one match is found, error out immediately
        if (count > 1) {
            print "ERROR: multiple matches found for file: " basename > "/dev/stderr"
            exit 1
        }
    }
    close(cmd)

    if (count == 0)
        return "NOT_FOUND"
    else
        return result
}

# Preserve header line (first row) unchanged
NR == 1 {
    print
    next
}

# Main processing block: runs once per data line
{
    b4 = basename($4)
    b5 = basename($5)

    new4 = find_unique(b4)
    new5 = find_unique(b5)

    # Log missing files to stderr
    if (new4 == "NOT_FOUND" || new5 == "NOT_FOUND") {
        print "WARNING (row " NR "): NOT_FOUND -> col4=" new4 ", col5=" new5 > "/dev/stderr"
    }

    $4 = new4
    $5 = new5

    print
}
' "$INPUT" > "$OUTPUT"

# Final confirmation message
echo "Done. Output written to $OUTPUT"

