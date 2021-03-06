#ifndef _XGC_DATA_READER_H
#define _XGC_DATA_READER_H

#include <string>
#include <fstream>
#include "io/xgcMesh.h"
#include "io/xgcData.h"
#include <adios.h>
#include <adios_read.h>
#include <adios_error.h>
#include <glob.h>
#include <mpi.h>

struct XGCDataReader {
  virtual ~XGCDataReader() {}

  virtual void open(const std::string&, bool filelist=false) = 0;
  virtual void close() = 0;
  virtual void read(XGCMesh& m, XGCData& d) = 0;
  virtual int advanceTimestep() {
    return ++ currentTimestep;
  }
  virtual int recedeTimestep() {
    return -- currentTimestep;
  }
 
  static std::vector<std::string> parseFileList(const std::string& filename) {
    std::vector<std::string> filenames;
    std::ifstream ifs(filename);
    std::string line;
    while (std::getline(ifs, line)) 
      filenames.push_back(line);
    ifs.close();
    return filenames;
  }

  static std::vector<std::string> glob(const std::string& pattern) {
    std::vector<std::string> filenames;
    glob_t results; 
    ::glob(pattern.c_str(), 0, NULL, &results); 
    for (int i=0; i<results.gl_pathc; i++)
      filenames.push_back(results.gl_pathv[i]); 
    globfree(&results);
    
    // fprintf(stderr, "find %lu files matching %s:\n", filenames.size(), pattern.c_str());
    // for (const auto &filename : filenames) 
    //   fprintf(stderr, "%s\n", filename.c_str());

    return filenames;
  }

  int getCurrentTimestep() const {return currentTimestep;}
  
  ADIOS_READ_METHOD adiosReadMethod() const {
    if (readMethod == "DATASPACES") return ADIOS_READ_METHOD_DATASPACES;
    else if (readMethod == "DIMES") return ADIOS_READ_METHOD_DIMES;
    else if (readMethod == "FLEXPATH") return ADIOS_READ_METHOD_FLEXPATH;
    else return ADIOS_READ_METHOD_BP;
  };

protected:
  MPI_Comm comm = MPI_COMM_WORLD;
  std::string readMethod;
  int currentTimestep = 0;

  friend class XGCDataReaderFactory;
};

struct XGCDataReader_Files : public XGCDataReader {
  std::vector<std::string> input_filename_list;

  void open(const std::string& pattern, bool filelist=false) {
    input_filename_list = XGCDataReader::glob(pattern);
    assert(input_filename_list.size() > 0);
  }

  void close() {}

  int advanceTimestep() {
    if (currentTimestep >= input_filename_list.size() - 1) return -1;
    else return ++ currentTimestep;
  }

  int recedeTimestep() {
    if (currentTimestep == 0) return -1;
    else return -- currentTimestep;
  }
};

struct XGCDataReader_AdiosStream : public XGCDataReader {
  ADIOS_FILE *varFP = NULL;

  void open(const std::string& name, bool filelist=false) {
    varFP = adios_read_open(name.c_str(), adiosReadMethod(), comm, ADIOS_LOCKMODE_ALL, -1.0);
  };

  void close() {
    // adios_close(*varFP);
  }

  void read(XGCMesh& m, XGCData& d) {
    d.readDpotFromADIOS(m, varFP);
  }

  int advanceTimestep() {
    adios_advance_step(varFP, 0, 1.0);
    return ++ currentTimestep;
  }

  int recedeTimestep() {return currentTimestep;}
};

struct XGCDataReader_AdiosFiles : public XGCDataReader_Files {
  void read(XGCMesh &m, XGCData& d) {
    ADIOS_FILE *varFP = adios_read_open_file(input_filename_list[currentTimestep].c_str(), adiosReadMethod(), comm);
    d.readDpotFromADIOS(m, varFP);
    // adios_close(*varFP);
  }
};

struct XGCDataReader_H5Files : public XGCDataReader_Files {
  void read(XGCMesh &m, XGCData& d) {
    d.readDpotFromH5(m, input_filename_list[currentTimestep]);
  }
};

struct XGCDataReaderFactory {
  static XGCDataReader* newXGCDataReader(const std::string& inputType, const std::string& readMethod, MPI_Comm comm) {
    XGCDataReader *r = NULL;
    if (inputType == "ADIOS_STREAM") r = new XGCDataReader_AdiosStream();
    else if (inputType == "ADIOS_FILES") r = new XGCDataReader_AdiosFiles();
    else if (inputType == "H5_FILES") r = new XGCDataReader_H5Files();
    else return r;

    r->comm = comm;
    r->readMethod = readMethod;
    return r;
  }
};

#endif
