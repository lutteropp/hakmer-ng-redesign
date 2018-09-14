/*
 * main.cpp
 *
 *  Created on: Sep 3, 2018
 *      Author: Sarah Lutteropp
 */

#include <iostream>

#include "external/CLI11.hpp"
#include "options.hpp"

int main(int argc, char* argv[]) {
	Options options;
	CLI::App app { "App description" };
	app.add_option("-f,--file", options.filepath, "A help string");
	CLI11_PARSE(app, argc, argv);
}
