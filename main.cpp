#include <iostream>
#include <functional>
#include <map>
#include <math.h>
#include <vector>
using mapstring = std::map<std::string,std::string>;
using functionstring = std::function<std::string(std::string,std::string)>;

int main(int argc, char **argv){
    using namespace std;
    map<string,function<string(vector<string>)>> map;
    map[""] = [](vector<string> w){
        return "Dodaj parametry:";
    };
    map["add"] = [](vector<string> w){
        string a = w.at(2);
        string b = w.at(3);
        a = to_string(stof(a)+stof(b));
        return a;};
    map["mod"] = [](vector<string> w){
        string a = w.at(2);
        string b = w.at(3);
        return to_string(stoi(a)%stoi(b));
    };
    map["sin"] = [](vector<string> w){
        string a = w.at(2);
        double sinus = stod(a);
        sinus = sin(sinus);
        return to_string(sinus);
    };
    vector<string> argumenty(argv, argv + argc);

    try{
        auto f = map.at(argumenty.at(1));
        cout<<f(argumenty);
    }catch (exception &e){
        cout<<"error\nDodaj parametry:";
    }
    return 0;
}