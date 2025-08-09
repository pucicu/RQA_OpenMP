/*************************************************************************/
/* RQA for very long data series                                         */
/* Norbert Marwan, Potsdam Institute for Climate Impact Research, 4/2009 */
/* License: GPLv3                                                        */
/*************************************************************************/

// compile for multithread with:  g++ -o rqa_omp rqa_omp.cpp -O3 -fopenmp


#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <iomanip>
#include <vector>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <omp.h>

//#define __VERBOSE__ // uncomment this line to get a more verbose output
    
enum norm { EUCLIDEAN, MAX, MIN, OP };
enum rptype { RP, CRP, JRP };
norm normType = EUCLIDEAN;
rptype rpType = RP;


/***********************************************************************************
 *  distance()  -  computes the Euclidean distance between 
 *                 phase space vectors X and Y which will be
 *                 reconstructed from a time-series by the
 *                 time delay method
 * 
 *         input:  pointer to the data vector U and pointer to 
 *                 the data vector V
 *         output: the Euclidean distance between the reconstructed
 *                 vectors X and Y.
 */
inline float dist(
    const std::vector<float>& data, 
    long i, 
    long j, 
    norm normType, 
    long m, 
    long t, 
    long cols
) {
        float d=0;
        float dNew;
        long k, k2;
        
        switch(normType) {
            case EUCLIDEAN:
                for ( k = 0; k < m; k++) {
                    for ( k2 = 0; k2 < cols; k2++) {
                        long idx_x = cols * (i + k * t) + k2;
                        long idx_y = cols * (j + k * t) + k2;
                        float diff = data[idx_x] - data[idx_y];
                        d += diff * diff;
                    }
                }
                d = sqrt(d);
                break;
            case MAX:
                d = 0;
                for ( k = 0; k < m; k++) {
                    for ( k2 = 0; k2 < cols; k2++) {
                        long idx_x = cols * (i + k * t) + k2;
                        long idx_y = cols * (j + k * t) + k2;
                        dNew = fabs(data[idx_x] - data[idx_y]);
                        if( dNew > d ) d = dNew;
                    }
                }
                break;
            case MIN:
                for ( k = 0; k < m; k++) {
                    for ( k2 = 0; k2 < cols; k2++) {
                        long idx_x = cols * (i + k * t) + k2;
                        long idx_y = cols * (j + k * t) + k2;
                        d += fabs(data[idx_x] - data[idx_y]);
                    }
                }
                break;
            case OP:
                d = fabs(data[i] - data[j]);
                break;
        }
        
        return d;
    }

/***********************************************************************************
 *  swap  -  swap two elements in a vector
 */
float swap(float& xPtr, float& yPtr)
{
    float temp;
    temp = xPtr;
    xPtr = yPtr;
    yPtr = temp;
    return 0;
}

/***********************************************************************************
 *  permute  -  permute the elements in a vector and assign the order with a symbol
 *              (create order pattern)
 */
long permute(long i1, float* x, long& symbol, long m, long cols)
{
    
    long output=0, output2=0, flag=0;
   
    for( long i = i1; i < m*cols; i++) {
        swap(x[i], x[i1]);
        if( i1 == m*cols - 1) {
            symbol++;
            flag = 1;
            // look for a specific order pattern
            for( long mi = 1; mi < m*cols; mi++) {
                if( x[mi - 1] <= x[mi] ) flag++; 
            }
            if( flag == m*cols ) {
                output = symbol;
            }
            
        } else {
            output2 = permute(i1+1, x, symbol, m, cols);
            if( output2 ) {
                output = output2;
                break;
            }
        }
        swap(x[i], x[i1]);
                
    }
    return output;
}

/***********************************************************************************
 *  makeOP  -  create a series of order patterns
 */
std::vector<float> makeOP(
    std::vector<float> datavector, 
    long m, 
    long t, 
    long cols
)
{
    // initialize variables and vectors
    long size = datavector.size();
    float *testVec;
    testVec = new float[m * cols];
    std::vector<float>::iterator dataPtr = datavector.begin();
    long symbol;
    #ifdef __VERBOSE__
    printf("Make Order Patterns.\n");
    #endif
    
    if(m * cols > 1) {
        
        for( long i = 0; i < size*cols ; i+=cols) {
            symbol  = 0;
            // make a vector of length m from the data series
            for( long mi = 0; mi < m; mi++) {
                for( long i_cols = 0; i_cols < cols; i_cols++) {
                    *(testVec + mi*cols + i_cols) = *(dataPtr + i + (mi * t) * cols + i_cols);
                }
            }
            // permute the elements and assign the order with the order pattern
            permute(0, testVec, symbol, m, cols);
            *(dataPtr + i) = float(symbol);
        }
        
    } else {
        std::cerr << "WARNING: Order patterns not created (dimension must be larger than one)." << std::endl;
    }

    return datavector;
}



/********************************************************************/
void usage(void)
{


  std::string version("$Revision: 1.1 $($Date: 2009/04/16 13:19:46 $)");
  
  int pos = 0;
  // remove "$"
  while ( pos >= 0  && pos < version.length()) {
      pos = version.find("$");
      if( pos >=0 && pos < version.length()) version.erase(pos, 1);
  }
  version.erase(0, 10); // remove "Revision:"
  pos = version.find("Date");
  version.erase(pos, 6); // remove "Date:"
  pos = version.find(" )");
  version.erase(pos, 1); // remove space before ")"

  std::cerr << "Commandline rqa multithreaded, version ";
  std::cerr << version << "\n";
  std::cerr << "  Calculates RR, DET and LMAX from a given file." << "\n";
  std::cerr << "  (c) Norbert Marwan, Potsdam Institute for Climate Impact Research (PIK)" << "\n\n";
  std::cerr << "usage:\n  rqa [options]\n\n";
  std::cerr << "options:\n";
  std::cerr << "    -i <string>   data filename (input)\n";
  std::cerr << "    -o <string>   filename RQA measures (output)\n";
  std::cerr << "    -p <string>   filename histogramme diagonal line lengths (output)\n";
  std::cerr << "    -q <string>   filename histogramme vertical line lengths (output)\n";
  std::cerr << "    -n <string>   distance norm (EUCLIDEAN, MAX, MIN, OP), default=EUCLIDEAN\n";
  std::cerr << "    -m <number>   embedding dimension, default=1\n";
  std::cerr << "    -t <number>   embedding delay, default=1\n";
  std::cerr << "    -e <number>   threshold, default=1\n";
  std::cerr << "    -l <number>   l_min, default=2\n";
  std::cerr << "    -v <number>   v_min, default=2\n";
  std::cerr << "    -w <number>   Theiler window, default=1\n";
  std::cerr << "    -s            silent (no messages displayed)\n";
  std::cerr << "    -h            print this help text\n\n";
}



/********************************************************************/
int main(int argc, char *argv[]) {
    
    // filenames and file-related info
    const char* inFilename = "none";
    const char* rqaFilename = "none";
    const char* delimiter_in = "\t";
    const char* delimiter = "\t";
    const char* histFilenameL = "none";
    const char* histFilenameV = "none";
    char ch, ch_old = *delimiter_in;
    std::ofstream fid;
    long cols = 0;
    float temp;
    
    // data vector
    std::vector<float> data;
    
    // embedding & RP parameters
    double e = 1;
    long m = 1;
    long t = 1;
    long tw=1;
    long lmin = 2;
    long vmin = 2;
    
    // RQA measures
    long RR = 0;
    double trueRR;
    double DET;
    double L;
    double LAM;
    double TT;
    long LMAX = 0;
    long VMAX = 0;
    double ENT = 0;
    double VENT = 0;
    double DIV = 0;

    // histograms of line lengths
    long *histv;
    long *histl;
    
    // time measurement
    time_t t1 = time(NULL);
    time_t t2;
    tm* systime;
    double elapsed_time;
    double elapsed_time_temp; // for each process separately
    
    // indices
    long k;
    long i; 
    long j;
    
    
    // helper variables for RQA
    long numberL = 0; // total number of diagonal lines
    long numberV = 0; // total number of vertical lines
    long line;
    long countL = 0; 
    long countV = 0; 
    double distance, d, oldDistance, distance2;
    long N;
    long cnt;
    long vRR = 0;
    
    // programme flags
    long silent = 0;
    long fileout = 0;

        
    /****************************************************************/
    // read the commandline input
    const char* optstring("i:o:p:q:n:m:t:e:l:w:v:sh?");
    std::string inString;

    if( argc==1 ) { usage(); return 0; }
    
    while(optind < argc)
    {
      optarg=NULL;
      int c;
      while( (c = getopt(argc, argv, optstring)) != -1)
      {
        switch(c)
        {
          case 'p':
            histFilenameL = argv[optind-1];
            break;
          case 'q':
            histFilenameV = argv[optind-1];
            break;
          case 'n':
            inString = argv[optind-1];
            if( !inString.compare(0, 3, std::string("EUC")) || !inString.compare(0, 3, std::string("euc")) ) normType = EUCLIDEAN;
            if( !inString.compare(0, 3, std::string("MAX")) || !inString.compare(0, 3, std::string("max")) ) normType = MAX;
            if( !inString.compare(0, 3, std::string("MIN")) || !inString.compare(0, 3, std::string("min")) ) normType = MIN;
            if( !inString.compare(0, 3, std::string("OP")) || !inString.compare(0, 3, std::string("op")) ) normType = OP;
            break;
          case 'i':
            inFilename = argv[optind-1];
            break;
          case 'o':
            rqaFilename = argv[optind-1];
            break;
          case 'e':
            e = atof(optarg);
            break;
          case 'm':
            m = atoi(optarg);
            break;
          case 't':
            t = atoi(optarg);
            break;
          case 'w':
            tw = atoi(optarg);
            break;
          case 'l':
            lmin = atoi(optarg);
            break;
          case 'v':
            vmin = atoi(optarg);
            break;
          case 's':
            silent = 1;
            break;
          case 'h':
            usage();
            exit(0);
          case '?':
            usage();
            exit(0);
        }
      }
    }

    
    if( strcmp(rqaFilename, "none") ) fileout = 1;
    
    /****************************************************************/
    // read file
    std::ifstream fid_in(inFilename); // open file
    
    // die if data file doesn't exist
    if( !fid_in ) {
        std::cerr << "ERROR: could not open file " << inFilename << std::endl;
        return 98;
    }
    


    std::string fline;
    bool first_line = true;
    while (std::getline(fid_in, fline)) {
        if (fline.empty()) continue; // skip blank lines
        std::istringstream iss(fline);
        float val;
        long col_count = 0;
        // Read values in this line
        while (iss >> val) {
            data.push_back(val);
            ++col_count;
        }
        if (first_line) {
            cols = col_count;
            first_line = false;
        } else if (col_count != cols) {
            std::cerr << "ERROR: Inconsistent column count in file at line: " << line << std::endl;
            return 99;
        }
    }
    fid_in.close();

    // length of the embedded time series                
    N = data.size()/cols - t * (m - 1);
    
    // make OP
    if( normType == OP ) {
        data = makeOP(data, m, t, cols);
        e = 0.001;
    }

    /****************************************************************/
    // verbose summary of input and settings
    if(!silent) {
#ifdef __PCWIN__
        if( rpType != RP ) {
            //std::cout << "Used files: " << inFilename << " and " << inFilename2 << " (" << i << " data points in " << cols << " column";
        } else {
            std::cout << "Used file: " << inFilename << " (" << data.size() << " data points in " << cols << " column";
        }
#else
        if( rpType != RP ) {
            //std::cout << "Used files: \033[1m" << inFilename << "\033[0m and \033[1m" << inFilename2 << "\033[0m (" << i << " data points in " << cols << " column";
        } else {
            std::cout << "Used file: \033[1m" << inFilename << "\033[0m (" << data.size() << " data points in " << cols << " column";
        }
#endif        
        if(cols > 1) std::cout << "s";
        std::cout << " read)." << std::endl;
        std::cout << "Used parameters: embedding dimension       m = " << m << std::endl;
        std::cout << "                 embedding delay           t = " << t << std::endl;
        std::cout << "                 recurrence threshold      e = " << e << std::endl;
        std::cout << "                 Theiler window            w = " << tw << std::endl;
        std::cout << "                 minimal diagonal line l_min = " << lmin << std::endl;
        std::cout << "                 minimal vertical line v_min = " << vmin << std::endl;
        std::cout << "                 ";
        switch(normType) {
            case EUCLIDEAN: std::cout << "Euclidean norm"; break;
            case MAX: std::cout << "maximum norm"; break;
            case MIN: std::cout << "minimum norm"; break;
            case OP: std::cout << "order patterns"; break;
        }
        std::cout << " used ";
        switch( rpType ) {
            case RP: std::cout << std::endl; break;
            case CRP: std::cout << "(compute cross recurrence plot)" << std::endl; break;
            case JRP: std::cout << "(compute joint recurrence plot)" << std::endl; break;
        }
    }
    
        
    /****************************************************************/
    // open results file
    if( fileout ) {
        fid.open(rqaFilename,std::ios::out);
        if( ! fid ) { 
            std::cerr << "ERROR: could not open " << rqaFilename << std::endl;
            fileout = 0;
        }
    }

    /****************************************************************/
    // RQA
                
    
    // initialise some variables and histograms
    RR = 0; cnt = 0;
    std::vector<float>::iterator dataPtr = data.begin();
    
    histl = new long[N+1];
    memset(histl, '\0', sizeof(long) * (N+1));
    if( histl == NULL ) {
        std::cerr << "Memory allocation for HISTL failed!" << std::endl;
        return 97;
    }
    histv = new long[N+1];
    memset(histv, '\0', sizeof(long) * (N+1));
    if( histv == NULL ) {
        std::cerr << "Memory allocation for HISTV failed!" << std::endl;
        return 95;
    }

            
    // separate the time series in segments for each process
    long i_start = 0;
    long i_end = N;
    
        
    /****************************************************************/
    LMAX = 0; 
    #pragma omp parallel if(N*m > 2000) 
    {
      #ifdef __VERBOSE__
        #ifdef _OPENMP
        if(omp_get_thread_num() == 0)
           std::cout <<  "                 Calculation is using \033[1m" << omp_get_num_threads() << "\033[0m processes." << std::endl;
        #endif
      #endif
    #pragma omp for private(i,j,k,d,distance,line,oldDistance) schedule(guided) reduction(+:cnt,RR)
       for( i = i_start; i < i_end; i++) {
           
           
            /****************************************************************/
            // diagonal-wise loop
            oldDistance = 0; line = 0; 
            // test for Theiler window
            if(i >= tw) {
                for( j = 0; j < N-i; j++) {

                    cnt++;

                    // calculate distance
                    distance = dist(data, i+j, j, normType, m, t, cols);
                    distance2 = distance;

                    // apply threshold
                    if (distance <= e) distance = 1.; else distance = 0.;

                    if(distance) 
//std::cout << "distance: " << distance2 << "\n";
//std::cout << "RP: " << distance << "\n";
//std::cout << "e: " << e << "\n";
                        RR++; // count all recurrences
                    if(oldDistance) line++; // count diagonal length of diagonal lines
                    else {
                        if(line) { 
                            #pragma omp atomic
                            histl[line]++; // make histogram of diagonal lines
                        }
                        line = 0;
                    }
                    oldDistance = distance;
                }
                if(line) 
                    #pragma omp atomic
                    histl[line]++; // make histogram of diagonal lines
            }
            /****************************************************************/
            // vertical-wise loop
            oldDistance = 0; line = 0; 
            for( j = 0; j < N; j++) {
                
                // calculate distance
                distance = dist(data, i, j, normType, m, t, cols);
        

                // apply threshold
                if (distance <= e) distance = 1.; else distance = 0.;

                if(oldDistance) line++; // count vertical length of vertical lines
                else {
                    if(line) { 
                        #pragma omp atomic
                        histv[line]++; // make histogram of vertical lines
                    }
                    line = 0;
                }
                oldDistance = distance;
            }
            if(line) 
                #pragma omp atomic
                histv[line]++; // make histogram of vertical lines
            /****************************************************************/
        }
    }
    
 


    /****************************************************************/
    
    
        
//         for(i=0; i<p; i++) {
//             std::cout << "\n LMAX: " << Lmax[i].l << " " << Lmax[i].pos << "\n";
//         }
        
        // recurrence rate
        trueRR = double(RR)/double(cnt);
        
        
        // count number of rec. points on diagonal lines
        // count the total number of diagonal lines
        // determine the longest diagonal line
        countL = 0; numberL = 0;
        for( i = lmin; i <= N-tw; i++) { 
            countL += histl[i] * (i);
            numberL += histl[i];
            if(histl[i]) LMAX = i;
        }
        
        // count number of rec. points on vertical lines
        // count the total number of vertical lines
        // determine the longest vertical line
        for( i = vmin; i <= N; i++) { 
            countV += histv[i] * (i);
            numberV += histv[i];
            if(histv[i]) VMAX = i;
        }
        
        // count all recurrence points for LAM
        for( i = 1; i <= VMAX; i++) vRR += histv[i] * i;
        
        // calculate probability
        double prob[N+1];  // prob. of a diagonal line with exact length L
        for( i = 1; i <= LMAX; i++) { 
            prob[i] = double(histl[i])/double(numberL);
        }
        
        // calculate L entropy
        for( i = lmin; i <= LMAX; i++) { 
            if(prob[i] > 0) ENT -= prob[i] * log(prob[i]);
        }
        
        // calculate probability
        for( i = 1; i <= VMAX; i++) { 
            prob[i] = double(histv[i])/double(numberV);
        }
        
        // calculate V entropy
        for( i = vmin; i <= VMAX; i++) { 
            if(prob[i] > 0) VENT -= prob[i] * log(prob[i]);
        }
        
        // calculate determinism and mean diagonal line length
        DET = double(countL) / double(RR);
        L = double(countL) / double(numberL);
        // calculate laminarity and trapping time
        LAM = double(countV) / double(vRR);
        TT = double(countV) / double(numberV);
        if(LMAX > 0) DIV = 1/double(LMAX); else DIV = 0;


        /****************************************************************/
        // output results into a file or on screen

        if(!silent) { 
            /****************************************************************/
            t2 = time(NULL) - t1 ;
            systime = localtime(&t2);
            std::cout.precision(6);
            std::cout << "\nComputation time: ";
            if( (systime->tm_mday - 1) ) std::cout << systime->tm_mday - 1 <<" day ";
            if( (systime->tm_hour - 1) > 0 ) std::cout << systime->tm_hour - 1 <<" h ";
            std::cout << systime->tm_min <<" min "<< systime->tm_sec<< " sec\n" << std::endl;
        }

        if( fileout ) {
        #ifdef __PCWIN__
                if(!silent) std::cout << "Write RQA results to file: " << rqaFilename ;
        #else
                if(!silent) std::cout << "Write RQA results to file: \033[1m" << rqaFilename << "\033[0m" ;
        #endif
            fid << trueRR << delimiter << DET << delimiter << LMAX << delimiter;
            fid << L << delimiter << ENT << delimiter << LAM << delimiter << TT << delimiter << VMAX;
            fid << std::endl;
        } else {
            if(!silent) {
                std::cout.precision(4);
                std::cout.width(6);
                std::cout << "Recurrence quantification analysis:\n";
                std::cout << "      RR:     " << trueRR << "\n";
                std::cout << "      DET:    " << DET;
                std::cout << "       \tLAM:     " <<LAM << "\n";
                if(RR > 0) 
                std::cout << "      DET/RR: " << DET/trueRR; 
                else 
                std::cout << "      DET/RR: " << "NaN";
                if(DET > 0) 
                std::cout << "      \tLAM/DET: " << LAM/DET;
                else 
                std::cout << "      \tLAM/DET: NaN";
                std::cout << "\n";
                std::cout << "      L_max:  " << LMAX;
                std::cout << "           \tV_max:   " << VMAX;
                std::cout << "\n";
                std::cout << "      L_mean: " << L;
                std::cout << "        \tTT:      " << TT;
                std::cout << "\n";
                std::cout << "      L_entr: " << ENT;
                std::cout << "        \tV_entr:  " << VENT;
                std::cout << "\n";
                if(DIV > 0) 
                std::cout << "      DIV:    " << DIV; 
                else 
                std::cout << "      DIV:    " << "NaN";
                std::cout << std::endl;
                
            }
        }
    
    
        if(!silent) { std::cout << std::endl; }
        
        // close results file
        if( fileout ) {
            fid.close();
        }    
        
        
        // store the histogramme of the diagonal line lengths
        if( strcmp(histFilenameL, "none") ) {
            std::ofstream fid(histFilenameL);
            if( ! fid ) { 
                std::cerr << "ERROR: could not open " << histFilenameL << std::endl;
            } else {
                for( i = 1; i <= N-tw; i++) { 
                    fid << i << ' ' << histl[i] << '\n';
                }
            }
                fid.close();
        }

        // store the histogramme of the vertical line lengths
        if( strcmp(histFilenameV, "none") ) {
            std::ofstream fid(histFilenameV);
            if( ! fid ) { 
                std::cerr << "ERROR: could not open " << histFilenameV << std::endl;
            } else {
                for( i = 1; i <= N; i++) { 
                    fid << i << ' ' << histv[i] << '\n';
                }
                fid.close();
            }
        }

    
    return 1;
}
