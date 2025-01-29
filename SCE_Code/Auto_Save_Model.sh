#!/bin/bash

# List of file names to be copied and stored
FILES=("main.cpp" "Bacteria.cpp" "Bacteria.hpp" "ConfigParser.cpp" "ConfigParser.h" "Diffusion2D.cpp" "Diffusion2D.hpp" "driver.cpp" "driver.h" "Fungi.cpp" "Fungi.hpp" "Grid.cpp" "Grid.hpp" "hyphaeSegment.cpp" "hyphaeSegment.h" "Nodes.cpp" "Nodes.hpp" "TissueBacteria.cpp" "TissueBacteria.hpp" "TissueGrid.cpp" "TissueGrid.hpp" "resources/bacteria_M.cfg")
DEST_FOLDER="Model_Code"
ZIP_NAME="Model_Code.zip"
FINAL_DESTINATION="$1"

# Check if destination is provided
if [[ -z "$FINAL_DESTINATION" ]]; then #-z tests if the string has zero length [[]] is a test command
    echo "Please provide the final destination as an argument."
    exit 1
fi

# Create folder
mkdir $DEST_FOLDER

# Copy files to the new folder
for FILE in "${FILES[@]}" #[@] when paired with double quotes treats all elements in array as separate
do
   cp $FILE $DEST_FOLDER/
done

# Zip the folder
zip -r $ZIP_NAME $DEST_FOLDER

# Move the zipped folder to the final destination
mv $ZIP_NAME $FINAL_DESTINATION

# Print the location of where the zipped folder is stored
echo "The zipped folder has been stored at: $FINAL_DESTINATION/$ZIP_NAME"

# Remove the original folder
rm -r $DEST_FOLDER
