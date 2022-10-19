#include <iostream>
#include <functional>
#include <math.h>
#include <random>
#include <time.h>
using namespace std;
using domain_t = std::vector<double>;
std::random_device rd;
std::mt19937 mt_generator(rd());

double optimise(auto function, auto domain, int maxIterations=1000){
    clock_t start, end;
    double cpu_time_used;
    start = clock();
    std::uniform_real_distribution<double> dist(domain.at(0), domain.at(1));
    double lowest=function(domain.at(0), domain.at(1));
    for(int i=0;i<maxIterations-1;i++){
        double rand1=dist(mt_generator);
        double rand2=dist(mt_generator);
        double temp;
        temp=function(rand1,rand2);
        if(temp<lowest)lowest=temp;
    }
    
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    cout<<"Iterations: "<<maxIterations<<" Time needed: "<<cpu_time_used<<" Result: ";
    return lowest;
}

int main(){
    double result=0;
    
    auto beal_f = [](double x, double y) {
        double firstPart = pow((1.5-x+(x*y)),2);
        double secondPart = pow(2.25-x+(x*pow(y,2)),2);
        double thirdPart = pow(2.625-x+x*pow(y,3),2);
        return firstPart+secondPart+thirdPart;
    };
    auto himmel_f = [](double x, double y) { return pow(pow(x,2)+y-11,2) + pow(x+pow(y,2)-7,2); };
    auto tcamel_f = [](double x, double y) {
        double firstPart = 2*pow(x,2);
        double secondPart = 1.05*pow(x,4);
        double thirdPart = (pow(x,6))/6;
        double fourthPart = x*y;
        double fifthPart = pow(y,2);
        return (firstPart-secondPart+thirdPart+fourthPart+fifthPart);
    };
    vector<double> domain={-4.5,4.5};
    for(int i=0;i<20;i++){
        cout<<optimise(beal_f,domain,10000)<<endl;
    }
    return 0;
}