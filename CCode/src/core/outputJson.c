/*
Output data into a json file
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "core/parson/parson.h"
#include "core/outputJson.h"


void outputJsonData(
    char * outputPath,
    char * title,
    struct DataSetForJson * dataSets,
    int numberOfDataSets)
{
    JSON_Value* root_value = json_value_init_object();
    JSON_Object* root_object = json_value_get_object(root_value);

    //set plot title
    json_object_set_string(root_object, "title", title);

    //set plot data
    JSON_Value* dataSets_Value = json_value_init_object();
    JSON_Object* dataSets_Object = json_value_get_object(dataSets_Value);
    int i, j;
    int contains_nan = 0;
    for (j = 0; j < numberOfDataSets; ++j) {
        contains_nan = 0;
        JSON_Value* arrayForData_Value = json_value_init_array();
        JSON_Array* arrayForData_Array = json_value_get_array(arrayForData_Value);
        for (i = 0; i < dataSets[j].length; ++i){
            if (isnan(dataSets[j].data[i]) || isinf(dataSets[j].data[i])){
                contains_nan = 1;
            }
        }
        if (contains_nan == 0){
            for (i = 0; i < dataSets[j].length; ++i) {
                json_array_append_number(arrayForData_Array, dataSets[j].data[i]);
            }
            json_object_dotset_value(dataSets_Object, dataSets[j].name, arrayForData_Value);
        }
    }

    json_object_dotset_value(root_object, "data arrays", dataSets_Value);

    //output to file
    FILE* outputFile = fopen(outputPath, "w");
    char * serialized_string = json_serialize_to_string_pretty(root_value);

    fprintf(outputFile, "%s", serialized_string);

}



void outputJsonData_and_Parameters(
    char * outputPath,
    char * title,
    struct DataSetForJson * dataSets,
    int numberOfDataSets,
    struct ParameterList * parameterList,
    int numberOfParameters)
{
    JSON_Value* root_value = json_value_init_object();
    JSON_Object* root_object = json_value_get_object(root_value);

    //set plot title
    json_object_set_string(root_object, "title", title);

    //set plot data
    JSON_Value* dataSets_Value = json_value_init_object();
    JSON_Object* dataSets_Object = json_value_get_object(dataSets_Value);
    int i, j;
    for (j = 0; j < numberOfDataSets; ++j) {
        JSON_Value* arrayForData_Value = json_value_init_array();
        JSON_Array* arrayForData_Array = json_value_get_array(arrayForData_Value);
        for (i = 0; i < dataSets[j].length; ++i) {
            json_array_append_number(arrayForData_Array, dataSets[j].data[i]);
        }
        json_object_dotset_value(dataSets_Object, dataSets[j].name, arrayForData_Value);
    }
    json_object_dotset_value(root_object, "data arrays", dataSets_Value);


    //set parameter data
    JSON_Value* parameters_Value = json_value_init_object();
    JSON_Object* parameters_Object = json_value_get_object(parameters_Value);
    for (j = 0; j < numberOfParameters; ++j) {
        json_object_set_number(parameters_Object, parameterList[j].name, parameterList[j].value);
    }
    json_object_dotset_value(root_object, "parameters", parameters_Value);


    

    //output to file
    FILE* outputFile = fopen(outputPath, "w");
    char * serialized_string = json_serialize_to_string_pretty(root_value);

    fprintf(outputFile, "%s", serialized_string);

}


