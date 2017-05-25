/*
outputJson.h
*/


struct DataSetForJson{
    char * name;
    double * data;
    int length;
};

struct ParameterList{
    char * name;
    double value;
};

void outputJsonData(
    char * outputPath,
    char * title,
    struct DataSetForJson * dataSets,
    int numberOfDataSets);

void outputJsonData_and_Parameters(
    char * outputPath,
    char * title,
    struct DataSetForJson * dataSets,
    int numberOfDataSets,
    struct ParameterList * parameterList,
    int numberOfParameters);