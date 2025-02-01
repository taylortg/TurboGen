#ifndef BLADE3D
#define BLADE3D

#include <array>
#include <string>

class Blade3D {
   public:
    Blade3D();
    ~Blade3D() = default;
    Blade3D(const Blade3D& other);
    Blade3D& operator=(const Blade3D& other);

    std::array<std::string, 9> meridionalChannelNames = {
        "x_leading", "y_leading", "z_leading", "x_trailing", "z_trailing", "x_hub", "z_hub", "x_shroud", "z_shroud"};
};

#endif  // BLADE3D