# Datatwin - Mining Simulation

Here we can find the code for the Astay - Mining Simulation, all of this functions are going to be implemented based on classes, functional tests and templates (if needed).

# Dependencies

- c++14
- gcc 
- [CMake](https://cmake.org/download/) - 3.28.1
- [Boost](https://www.boost.org/doc/libs/1_84_0/tools/build/doc/html/index.html#bbv2.installation) - 1.81.0
- [rapidcsv](https://github.com/d99kris/rapidcsv) - 8.80
- [SQLiteCpp](https://github.com/SRombauts/SQLiteCpp) - 3.3.1
- [sqlite3](https://www.linuxcapable.com/install-sqlite-on-fedora-linux/) - 3.42.0
- [gtest](https://gist.github.com/Cartexius/4c437c084d6e388288201aadf9c8cdd5) - 1.13.0

### Install CMake

Install cmake dependencies.
```
sudo dnf install gcc gcc-c++ openssl-devel bzip2-devel libffi-devel zlib-devel wget make -y
```

Download the latest `cmake-x.tar.gz` from the provided URL.
```
cd ~/Downloads/cmake-3.28.1
./bootstrap
make
sudo make install
```

### Install boost 

boost can be directly downloaded from fedora repo.
```
sudo dnf install boost-devel
```

### Install SQLiteCpp

To install SQLiteCpp let's create the thirdparty directory and download the latest version of it.
```
mkdir thirdparty && cd thirdparty
git clone https://github.com/SRombauts/SQLiteCpp.git
```

### Install sqlite3

Sqlite3 can be directly downloaded from fedora repo.
```
sudo dnf install sqlite3
sudo dnf install libsqlite3x-devel
```
The second command is needed to download the `sqlite3.h` file.

### Install gtest

Gtest can be directly downloaded from fedora repo.
```
sudo dnf install gtest
sudo dnf install gtest-devel
```
The second command is needed to download the `sqlite3.h` file.

# Compilation

Let's follow the next pipeline:

---
**NOTE**
It's recommended to compile whit -j2 as argument, however you are able to define x up to your number of cores

---

```
git clone https://<user>@bitbucket.org/astaysystems/astay-mining-sim.git
cd astay_mining_sim
mkdir build && cd build
cmake ..
make -j<x>
```

# Testing
Select which test to run before to compile. Modify this [CMakeLists.txt](tests/CMakeLists.txt) Be careful, you can select as many tests you want, but you must leave the `test_init.cpp` file as uncomment.
Let's run the next commands:
```
cd ../bin
./astay_mining_sim_test
```