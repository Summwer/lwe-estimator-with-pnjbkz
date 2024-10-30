#include <stdio.h>
#include <zipconf.h>
#include <zip.h>
#include <stdlib.h>

#define ZIPFILE_NAME  "testzip.zip"
#define SRCFILE_NAME   "test.txt"
#define FILENAME_IN_ZIP "test.txt"

int
main(int argc, char *argv[])
{
	int iErr = 0;
	struct zip * zipfile = NULL;
	struct zip_source * srcfile = NULL;

	zipfile = zip_open(ZIPFILE_NAME, ZIP_CREATE | ZIP_TRUNCATE, &iErr);
	if(!zipfile)
	{
		printf("zip open failed:%d\n", iErr);

		exit(EXIT_FAILURE);
	}

	//open a file for zip source
	srcfile = zip_source_file(zipfile, SRCFILE_NAME, 0, -1);
	if(!srcfile)
	{
		printf(" open source file failed\n");
		zip_close(zipfile);
		exit(EXIT_FAILURE);
	}

	//add file
	zip_file_add(zipfile, FILENAME_IN_ZIP, srcfile, ZIP_FL_OVERWRITE);

	zip_close(zipfile);
	
	return 0;
}