# Parallel Password Cracker

A multi-threaded implementation of a password cracker that can find the password associated with an MD5 hash.

### Prerequisites

A C++ compiler that supports c++11

### Compiling

if you use the "mariasfirstimage" image on a WINLAB sandbox1 node you will have to source the openssl lib to the script for compiling.
just run the make_macos.sh script by typing this:

```
./make_macos.sh
```

## Running

The release executable is located at ./brute_forcer  
The command line arguments give two options for arguments, the first being an MD5 hash, and the second thread count.\
For example:
```
./brute_forcer a6fe881cecd3fb7660083aea35cce430 8
```

## Authors

* **Shounak Rangwala**
* **Imad-uddin Siddiqui**
* **Amod Deo**



