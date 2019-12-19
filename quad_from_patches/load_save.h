#ifndef LOAD_SAVE_H
#define LOAD_SAVE_H

#include <string>
#include <vector>

std::vector<std::vector<size_t>> loadPatches(const std::string& filename);
std::vector<std::vector<size_t>> loadCorners(const std::string& filename);

#endif // LOAD_SAVE_H
