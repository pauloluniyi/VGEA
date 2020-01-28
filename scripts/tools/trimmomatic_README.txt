Trimmomatic can be run by either
(1) running the command
java -jar /path/to/my/trimmomatic/trimmomatic-XXX.jar
with the path and version number (XXX) replaced appropriately, or
(2) putting the shiver/tools/trimmomatic file in /path/to/my/trimmomatic/, then adding /path/to/my/trimmomatic/ to your PATH variable e.g. with
echo 'PATH=$PATH:/path/to/my/trimmomatic/' >> ~/.bashrc; source ~/.bashrc
then the simple command 'trimmomatic' will run trimmomatic, equivalent to the java -jar command above.
