#!/bin/bash
export MPP_DATA_DIRECTORY="$(cd "$(dirname "$0")/libs/mutationpp/data" && pwd)"
echo "MPP_DATA_DIRECTORY set to: $MPP_DATA_DIRECTORY"
./blast "$@"
