#include <algorithm>
#include <cassert>
#include "BVH.hpp"

BVHAccel::BVHAccel(std::vector<Object*> p, int maxPrimsInNode,
                   SplitMethod splitMethod)
    : maxPrimsInNode(std::min(255, maxPrimsInNode)), splitMethod(splitMethod),
      primitives(std::move(p))
{
    time_t start, stop;
    time(&start);
    if (primitives.empty())
        return;

    root = recursiveBuild(primitives);

    time(&stop);
    double diff = difftime(stop, start);
    int hrs = (int)diff / 3600;
    int mins = ((int)diff / 60) - (hrs * 60);
    int secs = (int)diff - (hrs * 3600) - (mins * 60);

    printf(
        "\rBVH Generation complete: \nTime Taken: %i hrs, %i mins, %i secs\n\n",
        hrs, mins, secs);
}

BVHBuildNode* BVHAccel::recursiveBuild(std::vector<Object*> objects)
{
    BVHBuildNode* node = new BVHBuildNode();

    // Compute bounds of all primitives in BVH node
    Bounds3 bounds;
    for (int i = 0; i < objects.size(); ++i)
        bounds = Union(bounds, objects[i]->getBounds());
    if (objects.size() == 1) {
        // Create leaf _BVHBuildNode_
        node->bounds = objects[0]->getBounds();
        node->object = objects[0];
        node->left = nullptr;
        node->right = nullptr;
        return node;
    }
    else if (objects.size() == 2) {
        node->left = recursiveBuild(std::vector{objects[0]});
        node->right = recursiveBuild(std::vector{objects[1]});

        node->bounds = Union(node->left->bounds, node->right->bounds);
        return node;
    }
    else {
        Bounds3 centroidBounds;
        for (int i = 0; i < objects.size(); ++i)
            centroidBounds =
                Union(centroidBounds, objects[i]->getBounds().Centroid());
        int dim = centroidBounds.maxExtent();
        switch (dim) {
        case 0:
            std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                return f1->getBounds().Centroid().x <
                       f2->getBounds().Centroid().x;
            });
            break;
        case 1:
            std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                return f1->getBounds().Centroid().y <
                       f2->getBounds().Centroid().y;
            });
            break;
        case 2:
            std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                return f1->getBounds().Centroid().z <
                       f2->getBounds().Centroid().z;
            });
            break;
        }
        if (splitMethod == SplitMethod::NAIVE)
        {
            auto beginning = objects.begin();
            auto middling = objects.begin() + (objects.size() / 2);
            auto ending = objects.end();

            auto leftshapes = std::vector<Object*>(beginning, middling);
            auto rightshapes = std::vector<Object*>(middling, ending);

            assert(objects.size() == (leftshapes.size() + rightshapes.size()));

            node->left = recursiveBuild(leftshapes);
            node->right = recursiveBuild(rightshapes);

            node->bounds = Union(node->left->bounds, node->right->bounds);
        }
        else    //SAH
        {
            double min_cost = std::numeric_limits<double>::max();
            int idx = -1;
            auto beginning = objects.begin();
            auto ending = objects.end();
            for (int i = 0; i < objects.size(); i++)
            {
                auto leftshapes = std::vector<Object*>(beginning, beginning + i);
                auto rightshapes = std::vector<Object*>(beginning + i, ending);
                Bounds3 left_bounds, right_bounds;
                for (int j = 0; j < leftshapes.size(); j++)
                {
                    left_bounds = Union(left_bounds, leftshapes[j]->getBounds().Centroid());
                }
                for (int j = 0; j < rightshapes.size(); j++)
                {
                    right_bounds = Union(right_bounds, rightshapes[j]->getBounds().Centroid());
                }
                auto S_left = left_bounds.SurfaceArea();
                auto a = leftshapes.size();
                auto S_right = right_bounds.SurfaceArea();
                auto b = rightshapes.size();
                auto S = bounds.SurfaceArea();
                
                auto cost = S_left / S * a + S_right / S * b;
                if (cost < min_cost)
                {
                    min_cost = cost;
                    idx = i;
                }
            }
            auto middling = objects.begin() + idx;

            auto leftshapes = std::vector<Object*>(beginning, middling);
            auto rightshapes = std::vector<Object*>(middling, ending);

            assert(objects.size() == (leftshapes.size() + rightshapes.size()));

            node->left = recursiveBuild(leftshapes);
            node->right = recursiveBuild(rightshapes);

            node->bounds = Union(node->left->bounds, node->right->bounds);
        }
    }

    return node;
}

Intersection BVHAccel::Intersect(const Ray& ray) const
{
    Intersection isect;
    if (!root)
        return isect;
    isect = BVHAccel::getIntersection(root, ray);
    return isect;
}

Intersection BVHAccel::getIntersection(BVHBuildNode* node, const Ray& ray) const
{
    // TODO Traverse the BVH to find intersection
    Intersection inter;
    std::array<int, 3>dirisNeg = { int(ray.direction.x > 0),int(ray.direction.y > 0),int(ray.direction.z > 0) };
    if (node->bounds.IntersectP(ray, ray.direction_inv, dirisNeg))
    {
        if (node->left == nullptr && node->right == nullptr)
        {
            return node->object->getIntersection(ray);
        }
        Intersection childL = BVHAccel::getIntersection(node->left, ray);
        Intersection childR = BVHAccel::getIntersection(node->right, ray);
        return childL.distance < childR.distance ? childL : childR;
    }
    return inter;
}