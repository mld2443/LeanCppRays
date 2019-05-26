///////////////
// RayTracer //
///////////////

#include <iostream>
#include <math.h>
#include <list>
#include <optional>
#include <random>

#include <png.h>

#define DTOR(x) x * 3.1415926535897932/180.0

using decimal = double;

//////////////////////
// Global Variables //
//////////////////////
std::random_device device;
std::mt19937 gen(device());
std::normal_distribution<decimal> nDist(0.0, 1.0);
std::uniform_real_distribution<decimal> uDist(0.0, 1.0);

//////////////////////////
// Forward Declarations //
//////////////////////////
class Material;
class Lambertian;
class Metallic;
class Dielectric;

class Shape;
class Plane;
class Sphere;

class Scene;
class Camera;

///////////////////
// MARK: Structs //
///////////////////
struct Range {
    decimal lower, upper;
    
    bool contains(const decimal& value) const { return lower <= value && value < upper; }
};


struct Vec3 {
    decimal x, y, z;
    
    decimal length() const;
    Vec3 normalize() const;
    bool isZero() const { return x == 0.0 && y == 0.0 && z == 0.0; }
};

Vec3 randomVector() { return Vec3{nDist(gen), nDist(gen), nDist(gen)}; }
Vec3 randomUnitVector() { return randomVector().normalize(); }

Vec3 operator+(const Vec3& lhs, const Vec3& rhs) { return Vec3{lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z}; }
Vec3 operator-(const Vec3& lhs, const Vec3& rhs) { return Vec3{lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z}; }
Vec3 operator-(const Vec3& v) { return Vec3{-v.x, -v.y, -v.z}; }
Vec3 operator*(const Vec3& lhs, const decimal rhs) { return Vec3{lhs.x * rhs, lhs.y * rhs, lhs.z * rhs}; }
Vec3 operator*(const decimal lhs, const Vec3& rhs) { return Vec3{lhs * rhs.x, lhs * rhs.y, lhs * rhs.z}; }
std::ostream& operator<<(std::ostream& str, const Vec3& v) { return str << "<" << v.x << ", " << v.y << ", " << v.z << ">"; }

Vec3 operator*(const Vec3& lhs, const Vec3& rhs) { return Vec3{lhs.x * rhs.x, lhs.y * rhs.y, lhs.z * rhs.z}; }
decimal dot(const Vec3& lhs, const Vec3& rhs) { return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z; }
Vec3 cross(const Vec3& lhs, const Vec3& rhs) { return Vec3{lhs.y * rhs.z - lhs.z * rhs.y, lhs.z * rhs.x - lhs.x * rhs.z, lhs.x * rhs.y - lhs.y * rhs.x}; }
Vec3 reflect(const Vec3& incoming, const Vec3& normal) { return incoming - (normal * (2.0 * dot(incoming, normal))); }
Vec3 refract(const Vec3& incoming, const Vec3& normal, const decimal eta) {
    const decimal cosI = -dot(incoming, normal), sinT2 = eta * eta * (1.0 - cosI * cosI);
    
    if (sinT2 > 1.0)
        return Vec3{};
    
    const decimal cosT = sqrt(1.0 - sinT2);
    return (incoming * eta) + (normal * (eta * cosI - cosT));
}

decimal Vec3::length() const { return sqrt(dot(*this, *this)); }
Vec3 Vec3::normalize() const { return *this * (1.0/length()); }


struct Ray {
    Vec3 origin, direction;
    
    Ray(const Vec3& o, const Vec3& d): origin(o), direction(d.normalize()) {}
    
    Vec3 project(const decimal dist) const { return origin + (direction * dist); }
};


struct Color {
    decimal r, g, b;
    
    Color(const decimal r = 0.0, const decimal g = 0.0, const decimal b = 0.0): r(r), g(g), b(b) {}
    Color(const char* desc) {
        unsigned int rr = 0, gg = 0, bb = 0;
        sscanf(desc, "#%2x%2x%2x", &rr, &gg, &bb);
        r = ((decimal)rr) / 255.0, g = ((decimal)gg) / 255.0, b = ((decimal)bb) / 255.0;
    }
    
    Color transform(const std::function<decimal(decimal)>& t) const { return Color{t(r), t(g), t(b)}; }
};

Color operator+(const Color& lhs, const Color& rhs) { return Color{lhs.r + rhs.r, lhs.g + rhs.g, lhs.b + rhs.b}; }
Color operator-(const Color& lhs, const Color& rhs) { return Color{lhs.r - rhs.r, lhs.g - rhs.g, lhs.b - rhs.b}; }
Color operator*(const Color& lhs, const Color& rhs) { return Color{lhs.r * rhs.r, lhs.g * rhs.g, lhs.b * rhs.b}; }
Color operator/(const Color& lhs, const decimal rhs) { return Color{lhs.r / rhs, lhs.g / rhs, lhs.b / rhs}; }
Color operator*(const Color& lhs, const decimal rhs) { return Color{lhs.r * rhs, lhs.g * rhs, lhs.b * rhs}; }
Color operator*(const decimal lhs, const Color& rhs) { return Color{lhs * rhs.r, lhs * rhs.g, lhs * rhs.b}; }


struct Pixel {
    png_byte r, g, b;
    
    Pixel(const png_byte r = 0, const png_byte g = 0, const png_byte b = 0): r(r), g(g), b(b) {}
    Pixel(const Color& c): r(clamp(c.r)), g(clamp(c.g)), b(clamp(c.b)) {}
    
private:
    static png_byte clamp(const decimal f) { return (png_byte) std::min(std::max((int) (f * 255.0), 0), 255); }
};


struct Intersection {
    decimal distance;
    Vec3 point;
    Vec3 normal;
    Material *material;
};


void abort_(const char * s, ...) {
    va_list args;
    va_start(args, s);
    vfprintf(stderr, s, args);
    fprintf(stderr, "\n");
    va_end(args);
    abort();
}


/////////////////////
// MARK: Materials //
/////////////////////
class Material {
public:
    Material(Color c): m_color(c) {}
    
    virtual Ray interact(const Ray&, const Vec3&, const Vec3&, const decimal) const=0;
    
    const Color m_color;
};


class Lambertian : public Material {
public:
    Lambertian(Color c): Material(c) {}
    
    Ray interact(const Ray& incoming, const Vec3& collision, const Vec3& normal, const decimal) const override {
        Vec3 target = collision + normal + randomUnitVector() * 0.999;
        return Ray{collision, target - collision};
    }
};


class Metallic : public Material {
public:
    Metallic(const Color c, const decimal f): Material(c), m_fuzz(f) {}
    
    Ray interact(const Ray& incoming, const Vec3& collision, const Vec3& normal, const decimal) const override {
        Vec3 reflected = reflect(incoming.direction, normal);
        
        if (m_fuzz > 0.0) {
            Vec3 fuzziness = randomUnitVector() * m_fuzz;
            decimal product = dot(fuzziness, normal);
            
            if (product < 0.0)
                fuzziness = fuzziness - (2.0 * product * normal);
            
            reflected = reflected + fuzziness;
        }
        
        return Ray{collision, reflected};
    }
private:
    decimal m_fuzz;
};


class Dielectric : public Material {
public:
    Dielectric(const Color c, const decimal i): Material(c.transform(sqrtl)), m_refractionIndex(i) {}
    
    Ray interact(const Ray& incoming, const Vec3& collision, const Vec3& normal, const decimal sceneIndex) const override {
        const decimal entering = dot(incoming.direction, normal);
        Vec3 refracted;
        
        if (entering > 0.0)
            refracted = refract(incoming.direction, -normal, m_refractionIndex / sceneIndex);
        else
            refracted = refract(incoming.direction, normal, sceneIndex / m_refractionIndex);

        if (refracted.isZero() || uDist(gen) < schlickApproximation(abs(entering), sceneIndex))
            return Ray{collision, reflect(incoming.direction, normal)};
        
        return Ray{collision, refracted};
    }
private:
    decimal m_refractionIndex;
    
    decimal schlickApproximation(const decimal cosX, const decimal sceneIndex) const {
        decimal r0 = (sceneIndex - m_refractionIndex) / (sceneIndex + m_refractionIndex);
        r0 *= r0;
        const decimal x = 1.0 - cosX;
        return r0 + (1.0 - r0) * x * x * x * x * x;
    }
};


//////////////////
// MARK: Shapes //
//////////////////
class Shape {
public:
    Shape(Material* m, const Vec3& p): m_position(p), m_material(m) {}
    
    std::optional<Intersection> intersectRay(const Ray& ray, const Range& window) {
        std::optional<decimal> distance = computeNearestIntersection(ray, window);
        
        if (!distance)
            return std::nullopt;
        
        Vec3 point = ray.project(*distance), normal = computeNormalAt(point);
        
        if (dot(ray.direction, normal) >= 0.0)
            return std::nullopt;
        
        return Intersection{*distance, point, normal, m_material};
    }
    
protected:
    virtual std::optional<decimal> computeNearestIntersection(const Ray&, const Range&) const=0;
    virtual Vec3 computeNormalAt(const Vec3&) const=0;

    const Vec3 m_position;
    Material *m_material;
};


class Plane : public Shape {
public:
    Plane(Material* m, const Vec3& p, const Vec3& n): Shape(m, p), m_normal(n.normalize()), m_normDotPos(dot(m_normal, p)) {}
    
protected:
    std::optional<decimal> computeNearestIntersection(const Ray& ray, const Range& window) const override {
        const decimal denominator = dot(m_normal, ray.direction);
        
        if (denominator == 0.0)
            return std::nullopt;
        
        const decimal distance = (m_normDotPos - dot(m_normal, ray.origin)) / denominator;
        
        if (window.contains(distance))
            return distance;
        
        return std::nullopt;
    }
    
    Vec3 computeNormalAt(const Vec3&) const override {
        return m_normal;
    }
    
private:
    const Vec3 m_normal;
    const decimal m_normDotPos;
};


class Sphere : public Shape {
public:
    Sphere(Material* m, const Vec3& p, const decimal r): Shape(m, p), m_radius(r * r) {}
    
protected:
    std::optional<decimal> computeNearestIntersection(const Ray& ray, const Range& window) const override {
        const Vec3 rCam = ray.origin - m_position;
        const Vec3 rRay = ray.direction;
        const decimal A = dot(rRay, rRay);
        const decimal B = dot(rCam, rRay);
        const decimal C = dot(rCam, rCam) - m_radius;
        const decimal square = B * B - A * C;
        
        if (square < 0.0)
            return std::nullopt;
        
        const decimal root = sqrt(square);
        const decimal D1 = (-B - root) / A;
        const decimal D2 = (-B + root) / A;
        
        if (window.contains(D1))
            return D1;
        if (window.contains(D2))
            return D2;
        
        return std::nullopt;
    }
    
    Vec3 computeNormalAt(const Vec3& point) const override {
        return (point - m_position).normalize();
    }
    
private:
    const decimal m_radius;
};


/////////////////
// MARK: Scene //
/////////////////
class Scene {
public:
    Scene(const Color h = "#4C7FFF",
          const Color s = "#FFFFFF",
          const decimal r = 1.0):
    m_horizon(h),
    m_sky(s),
    m_refractionIndex(r) {}
    
    ~Scene() {
        for (Shape* s : m_things) {
            if (s) {
                free(s);
            }
        }
    }
    
    void addShape(Shape* s) {
        m_things.push_back(s);
    }
    
    Color castRay(const Ray& start, const Range& frustum, const unsigned int depth) const {
        Ray ray = start;
        Range window = frustum;
        std::list<Color> colors{};
        
        while (true) {
            if (colors.size() >= depth)
                return Color{};
            
            std::optional<Intersection> nearest = findNearest(ray, window);
            
            if (!nearest)
                break;
            
            colors.push_back(nearest->material->m_color);
            ray = nearest->material->interact(ray, nearest->point, nearest->normal, m_refractionIndex);
            window = {window.lower, window.upper - nearest->distance};
        }
        
        Color result = skyBox(ray.direction);
        
        for (Color c : colors)
            result = result * c;
        
        return result;
    }
    
private:
    const Color m_horizon;
    const Color m_sky;
    const decimal m_refractionIndex;
    std::list<Shape*> m_things;
    
    std::optional<Intersection> findNearest(const Ray& ray, const Range& window) const {
        std::optional<Intersection> nearest = std::nullopt;
        Range currentWindow = window;
        
        for (Shape* s : m_things) {
            if (nearest)
                currentWindow = {window.lower, nearest->distance};
            
            std::optional<Intersection> candidate = s->intersectRay(ray, currentWindow);
            
            if (candidate) {
                nearest = candidate;
            }
        }
        
        return nearest;
    }
    
    Color skyBox(const Vec3& direction) const {
        const decimal interpolate = (0.5 * (direction.z + 1.0));
        return (m_horizon * (1.0 - interpolate)) + (m_sky * interpolate);
    }
};


//////////////////
// MARK: Camera //
//////////////////
class Camera {
public:
    Camera(const Vec3& position,
           const unsigned int width, const unsigned int height,
           const unsigned int sampling, const unsigned int depth,
           const Vec3& direction,
           const decimal FOV,
           const Vec3& up = Vec3{0.0, 0.0, 1.0}):
    m_position(position),
    m_width(width), m_height(height),
    m_sampling(sampling), m_depth(depth),
    m_frustum(Range{0.1, 1000.0}) {
        const Vec3 unitDirection = direction.normalize();
        const decimal screenWidth = tan(DTOR(FOV / 2.0));
        const decimal screenHeight = (((decimal) height) / ((decimal) width)) * screenWidth;
        const Vec3 iStar = cross(up, unitDirection).normalize();
        const Vec3 jStar = cross(iStar, unitDirection).normalize();
        
        m_iHat = iStar * (2.0 * screenWidth / (decimal) width);
        m_jHat = jStar * (2.0 * screenHeight / (decimal) height);
        m_origin = unitDirection + (iStar * -screenWidth) + (jStar * -screenHeight);
        
        m_film = new Pixel*[height];
        for (unsigned int y = 0u; y < height; ++y)
            m_film[y] = new Pixel[width];
    }
    
    ~Camera() {
        for (unsigned int y = 0u; y < m_height; ++y)
            free(m_film[y]);
        free(m_film);
    }
    
    void captureScene(const Scene& scene) {
        for (unsigned int y = 0u; y < m_height; y++) {
            for (unsigned int x = 0u; x < m_width; x++)
                this->m_film[y][x] = this->getPixel(scene, x, y);
            
            std::cout << ".";
        }
        
        std::cout << std::endl;
    }
    
    void developFilm(FILE* file) const {
        if (!file)
            abort_("[write_png_file] File could not be opened for writing");
        
        png_byte** image = (png_byte**) m_film;
        
        png_struct* negative = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
        if (!negative)
            abort_("[write_png_file] png_create_write_struct failed");
        
        png_info* info = png_create_info_struct(negative);
        if (!info)
            abort_("[write_png_file] png_create_info_struct failed");
        
        if (setjmp(png_jmpbuf(negative)))
            abort_("[write_png_file] Error during init_io");
        png_init_io(negative, file);
        
        
        if (setjmp(png_jmpbuf(negative)))
            abort_("[write_png_file] Error during writing header");
        png_set_IHDR(negative, info, m_width, m_height, 8, 2, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
        png_write_info(negative, info);
        
        if (setjmp(png_jmpbuf(negative)))
            abort_("[write_png_file] Error during writing bytes");
        png_write_image(negative, image);
        
        if (setjmp(png_jmpbuf(negative)))
            abort_("[write_png_file] Error during end of write");
        png_write_end(negative, NULL);
    }

private:
    const Vec3 m_position;
    const unsigned int m_width, m_height, m_sampling, m_depth;
    const Range m_frustum;
    Vec3 m_iHat, m_jHat, m_origin;
    Pixel** m_film;
    
    Pixel getPixel(const Scene& scene, const unsigned int x, const unsigned int y) const {
        Color sample{};
        
        for (unsigned int s = 0u; s < m_sampling; ++s) {
            const decimal xCoord = x + uDist(gen);
            const decimal yCoord = y + uDist(gen);
            
            const Vec3 screenSpacePosition = m_origin + (m_iHat * xCoord) + (m_jHat * yCoord);
            const Ray cast{m_position, screenSpacePosition};
            
            sample = sample + scene.castRay(cast, m_frustum, m_depth);
        }
        
        sample = sample / m_sampling;
        
        // Post-processing, makes the result brighter
        //sample = sample.transform(sqrt);
        
        return Pixel(sample);
    }
};


/////////////////////////
// MARK: Program Entry //
/////////////////////////
int main(int argc, const char * argv[]) {
    // Lights
    Scene s{};                     // Daytime lighting
    //Scene s{"#AA00AA", "#100460"}; // Nighttime lighting
    
    // Make-up
    Material* white = new Lambertian{"#FFFFFF"};
    Material* glass = new Dielectric{"#F0FFF0", 1.37};
    Material* matteCyan = new Lambertian{"#00FFFF"};
    Material* shinyYellow = new Metallic{"#FFFF00", 0.0};
    Material* gunmetal = new Metallic{"#808080", 0.1};
    Material* matteMagenta = new Lambertian{"#FF00FF"};
    
    // Actors
    s.addShape(new Plane{white, Vec3{}, Vec3{0.0,0.0,1.0}});
    s.addShape(new Sphere{glass, Vec3{9.0,0.0,6.0}, 6.0});
    s.addShape(new Sphere{matteCyan, Vec3{19.0,-5.0,3.0}, 3.0});
    s.addShape(new Sphere{shinyYellow, Vec3{20.0,5.0,4.0}, 4.0});
    s.addShape(new Sphere{gunmetal, Vec3{30.0,-16.0,16.0}, 16.0});
    s.addShape(new Sphere{matteMagenta, Vec3{100.0,170.0,30.0}, 30.0});
    
    // Camera
    //const unsigned int width = 320, height = 180, sampling = 100, depth = 10;
    const unsigned int width = 1920, height = 1080, sampling = 5000, depth = 10;
    Camera c = Camera(Vec3{0.0,0.0,7.0}, width, height, sampling, depth, Vec3{18.0,0.0,-1.0}, 100.0);
    
    // Action
    clock_t cycles = clock(); {
        c.captureScene(s);
    } cycles = clock() - cycles;
    
    std::cout << "Capture took " << cycles << " clocks, " << ((decimal)cycles)/CLOCKS_PER_SEC << " seconds." << std::endl;
    
    // Cut!
    FILE *file = fopen("Output.png", "wb");
    c.developFilm(file);
    fclose(file);
    
    // That's a wrap
    free(white);
    free(glass);
    free(matteCyan);
    free(shinyYellow);
    free(gunmetal);
    free(matteMagenta);
    
    return 0;
}
