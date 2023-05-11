
#include "load_lwechal.h"

//Download files
size_t dl_req_reply(void *buffer, size_t size, size_t nmemb, void *user_p)
{
	FILE *fp = (FILE *)user_p;
	size_t return_size = fwrite(buffer, size, nmemb, fp);
	//cout << (char *)buffer << endl;
	return return_size;
}

//http GET request for file download.
CURLcode dl_curl_get_req(const std::string &url, std::string filename)
{

	const char* file_name = filename.c_str();
	char* pc = new char[1024];//long enough
	strcpy(pc, file_name);

	FILE *fp = fopen(pc, "wb");

	//curl initilization 
	CURL *curl = curl_easy_init();
	// curl return
	CURLcode res = curl_global_init(CURL_GLOBAL_ALL);
	if (curl)
	{
		//set header of curl
		struct curl_slist* header_list = NULL;
		header_list = curl_slist_append(header_list, "User-Agent: Mozilla/5.0 (Windows NT 10.0; WOW64; Trident/7.0; rv:11.0) like Gecko");
		curl_easy_setopt(curl, CURLOPT_HTTPHEADER, header_list);

		//0: not accept 1: accept
		curl_easy_setopt(curl, CURLOPT_HEADER, 0);

		//URL address
		curl_easy_setopt(curl, CURLOPT_URL, url.c_str());

		//set ssl
		curl_easy_setopt(curl, CURLOPT_SSL_VERIFYPEER, false);
		curl_easy_setopt(curl, CURLOPT_SSL_VERIFYHOST, false);

		//CURLOPT_VERBOSE: 1, debug detail
		curl_easy_setopt(curl, CURLOPT_VERBOSE, 0);

		curl_easy_setopt(curl, CURLOPT_READFUNCTION, NULL);

		//set data accept function
		curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, &dl_req_reply);
		curl_easy_setopt(curl, CURLOPT_WRITEDATA, fp);

		curl_easy_setopt(curl, CURLOPT_NOSIGNAL, 1);

		//set requst the longest duration
		// curl_easy_setopt(curl, CURLOPT_CONNECTTIMEOUT, 6); // set transport and time out time  
		curl_easy_setopt(curl, CURLOPT_TIMEOUT, 6);

		// open request
		res = curl_easy_perform(curl);
	}
	// free curl 
	curl_easy_cleanup(curl);
	// free source
	fclose(fp);

	return res;
}

bool isFileExists_ifstream(string name) {
    ifstream f(name.c_str());
    return f.good();
}


LWEchal* load_lwe_challenge(int n, double alpha_){
    /*
    Load LWE challenge from file or website.

    :param n: LWE dimension
    :param alpha: the *standard deviation* of the secret is alpha*q

    */
    int alpha = int(round(alpha_ * 1000));
    ostringstream os, fname;
    const char* path = "lwechallenge/";
    fname << path << n << "-" << setw(3) << setfill('0')<< alpha <<".txt";

    if(not isFileExists_ifstream(fname.str())){
        os << "https://www.latticechallenge.org/lwe_challenge/challenges/LWE_" << n << "_" << setw(3) << setfill('0')<< alpha <<".txt";
        string dl_get_url = os.str();

        if(NULL==opendir(path))
            mkdir(path,0775);
        dl_curl_get_req(dl_get_url, fname.str());

        cout<<"Download \""<<fname.str()<<"\" successfully!"<<endl;
    }
    else{
        cout<<"\""<<fname.str()<<"\" exists."<<endl;
    }

	
    ifstream  fin;
    fin.open(fname.str(),ios::in);
    
	
    LWEchal* lwechal = new LWEchal;

	fin>>lwechal->n;
	fin>>lwechal->m;
	fin>>lwechal->q;
	fin>>lwechal->alpha;

	
	
	vector<Z_NR<ZT>> c;
	char ch;
	fin >> ch;
	if(ch == '['){
		while (fin >> ch && ch != ']'){
			fin.putback(ch);
			c.resize(c.size() + 1);
			// print_vector(c,0,c.size());
			if (!(fin >> c.back())){
				c.pop_back();
				break;
			}
		}
	}
	// print_vector(c,0,c.size());
	lwechal->c = c;
	
	

	vector<vector<Z_NR<ZT>>> matrix;
	fin>>ch;
	if(ch == '['){
		int i = 0;
		matrix.resize(0);
		// cout<<ch<<endl;
		while(fin >> ch && ch == '['){
			matrix.resize(i+1);
			// cout<<i<<endl;
			while (fin >> ch && ch != ']')
			{
				fin.putback(ch);
				matrix[i].resize(matrix[i].size() + 1);
				// cout<<lwechal->A[i].size()<<endl;
				if (!(fin >> matrix[i].back()))
				{
					matrix[i].pop_back();
					break;
				}
				// cout<<lwechal->A[i][0]<<endl;
			}
			i++;
		}
	}
	// print_matrix(matrix);
	


	lwechal->A.resize(lwechal->m,lwechal->n);
	for(int i = 0; i < lwechal->m ; i++){
		for(int j = 0; j < lwechal->n; j++){
			lwechal->A[i][j] = matrix[i][j]; 
		}
	}

	// lwechal->A.print(os);
	// cout<<os.str()<<endl;
	fin.close();
    return lwechal;
}




