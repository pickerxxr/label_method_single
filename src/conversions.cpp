#include <string.h>

#include "conversions.hpp"

// Removes \n from the end of line
void FIXLINE(char *s)
{
    int len = (int)strlen(s) - 1;
    if (s[len] == '\n')
        s[len] = 0;
}

void convert_edgelist(std::string inputfile, Converter *converter)
{
    FILE *inf = fopen(inputfile.c_str(), "r");
    size_t bytesread = 0;
    size_t linenum = 0;
    if (inf == NULL) {
        std::cout << "Could not load:" << inputfile
        << ", error: " << strerror(errno) << std::endl;
    }

    std::cout << "Reading in edge list format!" << std::endl;
    char s[1024];
    while (fgets(s, 1024, inf) != NULL) {
        linenum++;
        if (linenum % 10000000 == 0) {
            std::cout << "Read " << linenum << " lines, "
            << bytesread / 1024 / 1024. << " MB" << std::endl;
        }
        FIXLINE(s);
        bytesread += strlen(s);
        if (s[0] == '#')
            continue; // Comment
            if (s[0] == '%')
                continue; // Comment

                char delims[] = "\t, ";
                char *t;
                t = strtok(s, delims);
                if (t == NULL) {
                    std::cout << "Input file is not in right format. "
                    << "Expecting \"<from>\t<to>\". "
                    << "Current line: \"" << s << "\"\n";
                }
                vid_t from = atoi(t);
                t = strtok(NULL, delims);
                if (t == NULL) {
                    std::cout << "Input file is not in right format. "
                    << "Expecting \"<from>\t<to>\". "
                    << "Current line: \"" << s << "\"\n";
                }
                vid_t to = atoi(t);

                if (from != to) {
                    converter->add_edge(from, to);
                }
    }
    fclose(inf);
}

void textconvert_edgelist(std::string inputfile, TextConverter *converter)
{
    FILE *inf = fopen(inputfile.c_str(), "r");
    size_t bytesread = 0;
    size_t linenum = 0;
    if (inf == NULL) {
        std::cout << "Could not load:" << inputfile
        << ", error: " << strerror(errno) << std::endl;
    }

    std::cout << "Reading in edge list format!" << std::endl;
    char s[1024];
    while (fgets(s, 1024, inf) != NULL) {
        linenum++;
        if (linenum % 10000000 == 0) {
            std::cout << "Read " << linenum << " lines, "
            << bytesread / 1024 / 1024. << " MB" << std::endl;
        }
        FIXLINE(s);
        bytesread += strlen(s);
        if (s[0] == '#')
            continue; // Comment
            if (s[0] == '%')
                continue; // Comment

                char delims[] = "\t, ";
                char *t;
                t = strtok(s, delims);
                if (t == NULL) {
                    std::cout << "Input file is not in right format. "
                    << "Expecting \"<from>\t<to>\". "
                    << "Current line: \"" << s << "\"\n";
                }
                vid_t from = atoi(t);
                t = strtok(NULL, delims);
                if (t == NULL) {
                    std::cout << "Input file is not in right format. "
                    << "Expecting \"<from>\t<to>\". "
                    << "Current line: \"" << s << "\"\n";
                }
                vid_t to = atoi(t);

                if (from != to) {
                    converter->add_edge(from, to);
                }
    }
    fclose(inf);
}



void convert(std::string basefilename, Converter *converter)
{
    std::cout << "converting `" << basefilename << "'";
    if (basefilename.empty())
        std::cout << "empty file name" << std::endl;
    if (converter->done()) {
        std::cout << "skip" << std::endl;
        return;
    }
    converter->init();

    convert_edgelist(basefilename, converter);

    converter->finalize();
}


