/***************************************************************************/
/* Copyright(C) 2021


The authors of

Reliable Feature-Line Driven Quad-Remeshing
Siggraph 2021


 All rights reserved.
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
****************************************************************************/

#ifndef LOAD_SAVE_H
#define LOAD_SAVE_H

#include <string>
#include <vector>
std::vector<std::vector<size_t>> loadPatches(const std::string& filename);
std::vector<std::vector<size_t>> loadCorners(const std::string& filename);
std::vector<std::pair<size_t,size_t>> LoadFeatures(const std::string &filename);
std::vector<size_t> loadFeatureCorners(const std::string& filename);

#endif // LOAD_SAVE_H
