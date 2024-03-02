#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <fortran_executable> <input_file>\n", argv[0]);
        return EXIT_FAILURE;
    }

    char command[1024];
    snprintf(command, sizeof(command), "%s < %s", argv[1], argv[2]);

    int result = system(command);
    if (result != 0) {
        fprintf(stderr, "Error: Command '%s' failed with return code %d.\n", command, result);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
