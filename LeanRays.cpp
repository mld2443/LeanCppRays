///////////////
// RayTracer //
///////////////

#include <iostream>
#include <math.h>
#include <list>
#include <random>

#include "png.h"

#define DTOR(x) x * 3.14159265f/180.f

//////////////////////
// Global Variables //
//////////////////////
std::random_device rd;
std::mt19937 gen(rd());
std::normal_distribution<float> nDist(0.f, 1.f);
std::uniform_real_distribution<float> uDist(0.f, 1.f);

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
struct Float3 {
    float x, y, z;
    
    float length() const;
    Float3 normalize() const;
    bool isZero() const { return x == 0.f && y == 0.f && z == 0.f; }
};

Float3 randomVector() { return Float3{nDist(gen), nDist(gen), nDist(gen)}; }
Float3 randomUnitVector() { return randomVector().normalize(); }

Float3 operator+(const Float3& lhs, const Float3& rhs) { return Float3{lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z}; }
Float3 operator-(const Float3& lhs, const Float3& rhs) { return Float3{lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z}; }
Float3 operator-(const Float3& v) { return Float3{-v.x, -v.y, -v.z}; }
Float3 operator*(const Float3& lhs, const float rhs) { return Float3{lhs.x * rhs, lhs.y * rhs, lhs.z * rhs}; }
Float3 operator*(const float lhs, const Float3& rhs) { return Float3{lhs * rhs.x, lhs * rhs.y, lhs * rhs.z}; }
std::ostream& operator<<(std::ostream& str, const Float3& v) { str << "<" << v.x << ", " << v.y << ", " << v.z << ">"; return str; }

Float3 operator*(const Float3& lhs, const Float3& rhs) { return Float3{lhs.x * rhs.x, lhs.y * rhs.y, lhs.z * rhs.z}; }
float dot(const Float3& lhs, const Float3& rhs) { return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z; }
Float3 cross(const Float3& lhs, const Float3& rhs) { return Float3{lhs.y * rhs.z - lhs.z * rhs.y, lhs.z * rhs.x - lhs.x * rhs.z, lhs.x * rhs.y - lhs.y * rhs.x}; }
Float3 reflect(const Float3& incoming, const Float3& normal) { return incoming - (normal * (2.f * dot(incoming, normal))); }
Float3 refract(const Float3& incoming, const Float3& normal, const float eta) {
    const float cosI = -dot(incoming, normal), sinT2 = eta * eta * (1.f - cosI * cosI);
    
    if (sinT2 > 1.f)
        return Float3{};
    
    const float cosT = sqrtf(1.f - sinT2);
    return (incoming * eta) + (normal * (eta * cosI - cosT));
}

float Float3::length() const { return sqrtf(dot(*this, *this)); }
Float3 Float3::normalize() const { return *this * (1.f/length()); }


struct Ray {
    Float3 origin, direction;
    
    Ray(const Float3& o, const Float3& d): origin(o), direction(d.normalize()) {}
    
    Float3 project(const float dist) const { return origin + (direction * dist); }
};


struct Color {
    float r, g, b;
    
    Color(const float r = 0.f, const float g = 0.f, const float b = 0.f): r(r), g(g), b(b) {}
    
    Color(const char* desc) {
        unsigned int rr = 0, gg = 0, bb = 0;
        sscanf(desc, "#%2x%2x%2x", &rr, &gg, &bb);
        r = ((float)rr) / 255.f, g = ((float)gg) / 255.f, b = ((float)bb) / 255.f;
    }
    
    Color transform(const std::function<float(float)>& t) const { return Color{t(r), t(g), t(b)}; }
};

Color operator+(const Color& lhs, const Color& rhs) { return Color{lhs.r + rhs.r, lhs.g + rhs.g, lhs.b + rhs.b}; }
Color operator-(const Color& lhs, const Color& rhs) { return Color{lhs.r - rhs.r, lhs.g - rhs.g, lhs.b - rhs.b}; }
Color operator*(const Color& lhs, const Color& rhs) { return Color{lhs.r * rhs.r, lhs.g * rhs.g, lhs.b * rhs.b}; }
Color operator/(const Color& lhs, const float rhs) { return Color{lhs.r / rhs, lhs.g / rhs, lhs.b / rhs}; }
Color operator*(const Color& lhs, const float rhs) { return Color{lhs.r * rhs, lhs.g * rhs, lhs.b * rhs}; }
Color operator*(const float lhs, const Color& rhs) { return Color{lhs * rhs.r, lhs * rhs.g, lhs * rhs.b}; }


png_byte channelToByte(const float f) { return (png_byte) std::min(std::max((int) (f * 255.0f), 0), 255); }

struct Pixel {
    png_byte r, g, b;
    
    Pixel(const png_byte r = 0, const png_byte g = 0, const png_byte b = 0): r(r), g(g), b(b) {}
    Pixel(const Color& c): r(channelToByte(c.r)), g(channelToByte(c.g)), b(channelToByte(c.b)) {}
};


struct Intersection {
    float distance;
    Float3 point;
    Float3 normal;
    Material *material;
};


template <typename T>
struct Range {
    T lower, upper;
    
    bool contains(const T& value) const { return lower <= value && value < upper; }
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
    
    virtual Ray interact(const Ray&, const Float3&, const Float3&, const float) const=0;
    
    const Color m_color;
};

class Lambertian : public Material {
public:
    Lambertian(Color c): Material(c) {}
    
    Ray interact(const Ray& incoming, const Float3& collision, const Float3& normal, const float) const {
        Float3 target = collision + normal + randomUnitVector() * 0.999f;
        return Ray{collision, target - collision};
    }
};

class Metallic : public Material {
public:
    Metallic(const Color c, const float f): Material(c), m_fuzz(f) {}
    
    Ray interact(const Ray& incoming, const Float3& collision, const Float3& normal, const float) const {
        Float3 reflected = reflect(incoming.direction, normal);
        
        if (m_fuzz > 0.f) {
            Float3 fuzziness = randomUnitVector() * m_fuzz;
            float product = dot(fuzziness, normal);
            
            if (product < 0.f)
                fuzziness = fuzziness - (2.f * product * normal);
            
            reflected = reflected + fuzziness;
        }
        
        return Ray{collision, reflected.normalize()};
    }
private:
    float m_fuzz;
};

class Dielectric : public Material {
public:
    Dielectric(const Color c, const float i): Material(c.transform(sqrtf)), m_refractionIndex(i) {}
    
    Ray interact(const Ray& incoming, const Float3& collision, const Float3& normal, const float sceneIndex) const {
        const float entering = dot(incoming.direction, normal);
        float cosX;
        Float3 refracted;
        
        if (entering > 0.f) {
            cosX = entering;
            refracted = refract(incoming.direction, -normal, m_refractionIndex / sceneIndex);
        } else {
            cosX = -entering;
            refracted = refract(incoming.direction, normal, sceneIndex / m_refractionIndex);
        }
        
        if (refracted.isZero() || uDist(gen) < schlickApproximation(cosX, sceneIndex))
            return Ray{collision, reflect(incoming.direction, normal)};
        
        return Ray{collision, refracted};
    }
private:
    float m_refractionIndex;
    
    float schlickApproximation(const float cosX, const float sceneIndex) const {
        float r0 = (sceneIndex - m_refractionIndex) / (sceneIndex + m_refractionIndex);
        r0 *= r0;
        const float x = 1.f - cosX;
        return r0 + (1.f - r0) * x * x * x * x * x;
    }
};

//////////////////
// MARK: Shapes //
//////////////////
class Shape {
public:
    Shape(Material* m, const Float3& p): m_position(p), m_material(m) {}
    
    virtual Float3 computeNormalAt(const Float3&) const=0;
    virtual float computeNearestIntersection(const Ray&, const Range<float>&) const=0;
    
    Intersection* intersectRay(const Ray& ray, const Range<float>& window) {
        float distance = computeNearestIntersection(ray, window);
        
        if (distance < 0.f)
            return nullptr;
        
        Float3 point = ray.project(distance), normal = computeNormalAt(point);
        
        if (dot(ray.direction, normal) >= 0.f)
            return nullptr;
        
        return new Intersection{distance, point, normal, m_material};
    }
    
protected:
    const Float3 m_position;
    Material *m_material;
};

class Plane : public Shape {
public:
    Plane(Material* m, const Float3& p, const Float3& n): Shape(m, p), m_normal(n.normalize()), m_normDotPos(dot(m_normal, p)) {}
    
    Float3 computeNormalAt(const Float3&) const {
        return m_normal;
    }
    
    float computeNearestIntersection(const Ray& ray, const Range<float>& window) const {
        const float denominator = dot(m_normal, ray.direction);
        
        if (denominator == 0.f)
            return -1.f;
        
        const float distance = (m_normDotPos - dot(m_normal, ray.origin)) / denominator;
        
        if (window.contains(distance))
            return distance;
        
        return -1.f;
    }
    
private:
    const Float3 m_normal;
    const float m_normDotPos;
};

class Sphere : public Shape {
public:
    Sphere(Material* m, const Float3& p, const float r): Shape(m, p), m_radius(r * r) {}
    
    Float3 computeNormalAt(const Float3& point) const {
        return (point - m_position).normalize();
    }
    
    float computeNearestIntersection(const Ray& ray, const Range<float>& window) const {
        const Float3 rCam = ray.origin - m_position;
        const Float3 rRay = ray.direction;
        const float A = dot(rRay, rRay);
        const float B = dot(rCam, rRay);
        const float C = dot(rCam, rCam) - m_radius;
        const float square = B * B - A * C;
        
        if (square < 0.f)
            return -1.f;
        
        const float root = sqrt(square);
        const float D1 = (-B - root) / A;
        const float D2 = (-B + root) / A;
        
        if (window.contains(D1))
            return D1;
        if (window.contains(D2))
            return D2;
        return -1.f;
    }
    
private:
    const float m_radius;
};

/////////////////
// MARK: Scene //
/////////////////
class Scene {
public:
    Scene(const Color h = {0.3f, 0.5f, 1.f},
          const Color s = {1.f, 1.f, 1.f},
          const float r = 1.f):
    horizon(h),
    sky(s),
    refractionIndex(r) {}
    
    ~Scene() {
        for (Shape* s : things)
            free(s);
    }
    
    void addShape(Shape* s) {
        things.push_back(s);
    }
    
    Color castRay(const Ray& ray, const Range<float>& frustum, const unsigned int depth) const {
        Ray dir = ray;
        Range<float> window = frustum;
        std::list<Color> colors{};
        
        while (true) {
            if (colors.size() >= depth)
                return Color{};
            
            Intersection* nearest = findNearest(dir, window);
            
            if (!nearest)
                break;
            
            colors.push_back(nearest->material->m_color);
            dir = nearest->material->interact(dir, nearest->point, nearest->normal, refractionIndex);
            window = {window.lower, window.upper - nearest->distance};
            
            free(nearest);
        }
        
        Color sky = skyBox(dir.direction);
        
        for (Color c : colors)
            sky = sky * c;
        
        return sky;
    }
    
private:
    const Color horizon;
    const Color sky;
    const float refractionIndex;
    std::list<Shape*> things;
    
    Intersection* findNearest(const Ray& ray, const Range<float>& window) const {
        Intersection* nearest = nullptr;
        Range<float> currentWindow = window;
        
        for (Shape* s : things) {
            if (nearest)
                currentWindow = {window.lower, nearest->distance};
            
            Intersection* candidate = s->intersectRay(ray, currentWindow);
            
            if (candidate) {
                if (nearest)
                    delete nearest;
                
                nearest = candidate;
            }
        }
        
        return nearest;
    }
    
    Color skyBox(const Float3& direction) const {
        const float interpolate = (0.5f * (direction.z + 1.f));
        return (horizon * (1.f - interpolate)) + (sky * interpolate);
    }
};

//////////////////
// MARK: Camera //
//////////////////
class Camera {
public:
    Camera(const Float3& position,
           const unsigned int width, const unsigned int height,
           const unsigned int sampling, const unsigned int depth,
           const Float3& direction,
           const float FOV,
           const Float3& up = Float3{0.f, 0.f, 1.f}):
    position(position),
    width(width), height(height),
    sampling(sampling), depth(depth),
    frustum({0.1f, 1000.f}) {
        const Float3 unitDirection = direction.normalize();
        const float screenWidth = tan(DTOR(FOV / 2.f));
        const float screenHeight = (((float) height) / ((float) width)) * screenWidth;
        const Float3 iStar = cross(up, unitDirection).normalize();
        const Float3 jStar = cross(iStar, unitDirection).normalize();
        
        iHat = iStar * (2.f * screenWidth / (float) width);
        jHat = jStar * (2.f * screenHeight / (float) height);
        origin = unitDirection + (iStar * -screenWidth) + (jStar * -screenHeight);
        
        film = new Pixel*[height];
        for (unsigned int y = 0; y < height; ++y)
            film[y] = new Pixel[width];
    }
    
    ~Camera() {
        for (unsigned int y = 0; y < height; ++y)
            free(film[y]);
        free(film);
    }
    
    void captureScene(const Scene& scene) {
        for (unsigned int y = 0; y < height; y++) {
            for (unsigned int x = 0; x < width; x++)
                film[y][x] = getPixel(scene, x, y);
            
            std::cout << ".";
        }
        
        std::cout << std::endl;
    }
    
    void developFilm(FILE* file) const {
        if (!file)
            abort_("[write_png_file] File could not be opened for writing");
        
        png_byte** image = (png_byte**) film;
        
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
        png_set_IHDR(negative, info, width, height, 8, 2, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
        png_write_info(negative, info);
        
        if (setjmp(png_jmpbuf(negative)))
            abort_("[write_png_file] Error during writing bytes");
        png_write_image(negative, image);
        
        if (setjmp(png_jmpbuf(negative)))
            abort_("[write_png_file] Error during end of write");
        png_write_end(negative, NULL);
    }

private:
    const Float3 position;
    const unsigned int width, height, sampling, depth;
    const Range<float> frustum;
    Float3 iHat, jHat, origin;
    Pixel** film;
    
    Pixel getPixel(const Scene& scene, const unsigned int x, const unsigned int y) const {
        Color sample{};
        
        for (unsigned int s = 0; s < sampling; ++s) {
            const float xCoord = x + uDist(gen);
            const float yCoord = y + uDist(gen);
            
            const Float3 screenSpacePosition = origin + (iHat * xCoord) + (jHat * yCoord);
            const Ray cast{position, screenSpacePosition};
            
            sample = sample + scene.castRay(cast, frustum, depth);
        }
        
        sample = sample / sampling;
        
        // Post-processing, makes the result brighter
        //sample = sample.transform(sqrtf);
        
        return Pixel(sample);
    }
};

/////////////////////////
// MARK: Program Entry //
/////////////////////////
int main(int argc, const char * argv[]) {
    // Lights
    Scene s{};                                    // Daytime lighting
    //Scene s{Color{"#AA00AA"}, Color{"#100460"}};  // Nighttime lighting
    
    Material* white = new Lambertian(Color{"#FFFFFF"});
    Material* glass = new Dielectric(Color{"#F0FFF0"}, 1.37f);
    Material* matteCyan = new Lambertian(Color{"#00FFFF"});
    Material* shinyYellow = new Metallic(Color{"#FFFF00"}, 0.f);
    Material* gunmetal = new Metallic(Color{"#808080"}, 0.1f);
    Material* matteMagenta = new Lambertian(Color{"#FF00FF"});
    
    s.addShape(new Plane(white, Float3{0.f,0.f,0.f}, Float3{0.f,0.f,1.f}));
    s.addShape(new Sphere(glass, Float3{9.f,0.f,6.f}, 6.f));
    s.addShape(new Sphere(matteCyan, Float3{19.f,-5.f,3}, 3.f));
    s.addShape(new Sphere(shinyYellow, Float3{20.f,5.f,4.f}, 4.f));
    s.addShape(new Sphere(gunmetal, Float3{30.f,-16.f,16.f}, 16.f));
    s.addShape(new Sphere(matteMagenta, Float3{100.f,170.f,30.f}, 30.f));
    
    // Camera
    const unsigned int width = 640, height = 360, sampling = 100, depth = 10;
    //const unsigned int width = 1280, height = 720, sampling = 1000, depth = 10;
    Camera c = Camera(Float3{0.f,0.f,7.f}, width, height, sampling, depth, Float3{18.f,0.f,-1.f}, 100.f);
    
    // Action
    clock_t cycles = clock(); {
        c.captureScene(s);
    } cycles = clock() - cycles;
    
    std::cout << "Capture took " << cycles << " clocks, " << ((float)cycles)/CLOCKS_PER_SEC << " seconds." << std::endl;
    
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
