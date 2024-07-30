#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
    if (argc != 2 && argc != 3) {
        fprintf(stderr, "Usage: %s <executable> <optional_input_file>\n", argv[0]);
        return EXIT_FAILURE;
    }

    char command[1024];
    if (argc == 2) {
        snprintf(command, sizeof(command), "%s", argv[1]);
    } else {
        snprintf(command, sizeof(command), "%s < %s", argv[1], argv[2]);
    }

    int result = system(command);
    if (result != EXIT_SUCCESS) {
        fprintf(stderr, "Error: Command '%s' failed with return code %d.\n", command, result);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
