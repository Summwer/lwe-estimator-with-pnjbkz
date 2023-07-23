
#include "utils.h"
#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <cstddef>
#include<string>
#include<curl/curl.h>

using namespace std;

//Download files
size_t dl_req_reply(void *buffer, size_t size, size_t nmemb, void *user_p);

//http GET request for file download.
CURLcode dl_curl_get_req(const std::string &url, std::string filename);

bool isFileExists_ifstream(string name);

LWEchal* load_lwe_challenge(int n, double alpha_);