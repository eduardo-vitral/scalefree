######################################################################
# Install the Scale-free modeling software
# Just give the command
#    prompt> source install.com
######################################################################

# Remove any existing executables and object files
\rm scalefree.o
\rm scalefree.e

# Compile the program
make all

# Write a message
echo ' '
echo 'Scale-free modeling software installed'
echo '  Documentation in file README'
echo '  Executable is scalefree.e'
echo '  Example input files in directory examp'
echo ' '
