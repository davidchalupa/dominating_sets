#include <stdio.h>
#include "cli.h"

int main()
{
    cli *cli_instance;

    cli_instance = new cli();
    cli_instance->start_cli();
    delete(cli_instance);

    return 0;
}
