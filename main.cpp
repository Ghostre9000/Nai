#include <any>
#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <cmath>
using mojamapa_t = std::map<std::string , std::string>;
using mojafunkcja_t = std::function<double(std::vector<double>)>;
void wypisz(mojamapa_t mapa, mojafunkcja_t fun) {
    using namespace std;
    for (auto kv : mapa) {
        auto [k, v] = kv;
        cout << "klucz: " << k << "; wartosc " << fun(v) << endl;
    }
}
int main(int argc, char **argv) {
    using namespace std;
    map<string , string> mapa = {{"Key", "Value"}};
    map<string, mojafunkcja_t> formatery;
    formatery["mod"] = [](vector<double> x) { return fmod(x[0],x[1]); };
    formatery["add"] = [](vector<double> x) { return x[0]+x[1]; };
    formatery["sin"] = [](vector<double> x) { return sin(x[0]); };
    try {
        vector<string> argumenty(argv, argv + argc);
        auto selected_f = argumenty.at(1);
        wypisz(mapa, formatery.at(selected_f));
    } catch (std::out_of_range aor) {
        cout << "podaj argument. Dostepne to: ";
        for (auto [k, v] : formatery) cout << " " << k;
        cout << endl;
    }
    return 0;
}