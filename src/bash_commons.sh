#!/bin/bash

# --- Global settings ----------------------------------------------------------

# Strict mode options
#set -e           # "exit-on-error" shell option
set -u           # "no-unset" shell option
#set -o pipefail  # exit on within-pipe error
#set -o errtrace  # ERR trap inherited by shell functions

# For a friendlier use of colors in Bash
red=$'\e[1;31m' # Red
grn=$'\e[1;32m' # Green
yel=$'\e[1;33m' # Yellow
blu=$'\e[1;34m' # Blue
mag=$'\e[1;35m' # Magenta
cya=$'\e[1;36m' # Cyan
end=$'\e[0m'

function _extract_mtpdb {
    # Local variables
    local db_filename="MTPDB.sqlite"
    local archive_path="./data/in/MTP-DB/${db_filename}.gz"
    local db_path="$1"

    # Extract the archive
    if [[ ! -f "$db_path" ]]; then
        if [[ -f "$archive_path" ]]; then
            printf "Extracting MTP-DB ... "
            gzip -dc "$archive_path" > "$db_path"
            printf "Done\n"
        else
            printf "\nMTP-DB not found locally!"
            printf "Run 'kerblam data fetch' to download fresh database.\n"
        fi
    fi
}
