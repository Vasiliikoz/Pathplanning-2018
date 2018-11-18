#include <iomanip>
#include <iostream>

using namespace std;

int main(int argc, char *argv[]) {
    if (argc == 1)
        return 0;
    auto per = argv[1];
    int64_t n = 0;
    while (static_cast<int>(*per) != 0) {
        n =(n * 10) + static_cast<int>(*per) - 48;
        per++;
    }
    while (n-- > 0)
        cout << "Pathplanning" << "\n";
}
