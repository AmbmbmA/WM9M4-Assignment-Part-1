#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>

#include "GamesEngineeringBase.h" // Include the GamesEngineeringBase header
#include <algorithm>
#include <chrono>

#include <cmath>
#include "matrix.h"
#include "colour.h"
#include "mesh.h"
#include "zbuffer.h"
#include "renderer.h"
#include "RNG.h"
#include "light.h"
#include "triangle.h"

#include <immintrin.h>
#include <thread>
#include <future>
// Utility function to generate a random rotation matrix
// No input variables
matrix makeRandomRotation() {
	RandomNumberGenerator& rng = RandomNumberGenerator::getInstance();
	unsigned int r = rng.getRandomInt(0, 3);

	switch (r) {
	case 0: return matrix::makeRotateX(rng.getRandomFloat(0.f, 2.0f * M_PI));
	case 1: return matrix::makeRotateY(rng.getRandomFloat(0.f, 2.0f * M_PI));
	case 2: return matrix::makeRotateZ(rng.getRandomFloat(0.f, 2.0f * M_PI));
	default: return matrix::makeIdentity();
	}
}

/*
void renderNew(Renderer& renderer, Mesh* mesh, matrix& camera, Light& L, matrix& p, size_t initindex, size_t endindex) {

	// Iterate through all triangles in the mesh
	for (size_t ind = initindex; ind < endindex; ind++) {
		Vertex t[3]; // Temporary array to store transformed triangle vertices

		// Transform each vertex of the triangle
		for (unsigned int i = 0; i < 3; i++) {
			t[i].p = p * mesh->vertices[mesh->triangles[ind].v[i]].p; // Apply transformations
			t[i].p.divideW(); // Perspective division to normalize coordinates

			// Transform normals into world space for accurate lighting
			// no need for perspective correction as no shearing or non-uniform scaling
			t[i].normal = mesh->world * mesh->vertices[mesh->triangles[ind].v[i]].normal;
			t[i].normal.normalise();

			// Map normalized device coordinates to screen space
			t[i].p[0] = (t[i].p[0] + 1.f) * 0.5f * static_cast<float>(renderer.canvas.getWidth());
			t[i].p[1] = (t[i].p[1] + 1.f) * 0.5f * static_cast<float>(renderer.canvas.getHeight());
			t[i].p[1] = renderer.canvas.getHeight() - t[i].p[1]; // Invert y-axis

			// Copy vertex colours
			t[i].rgb = mesh->vertices[mesh->triangles[ind].v[i]].rgb;
		}

		// Clip triangles with Z-values outside [-1, 1]
		if (fabs(t[0].p[2]) > 1.0f || fabs(t[1].p[2]) > 1.0f || fabs(t[2].p[2]) > 1.0f) continue;

		// Create a triangle object and render it
		triangle tri(t[0], t[1], t[2]);
		tri.draw(renderer, L, mesh->ka, mesh->kd);
	}
}

void renderMT(Renderer& renderer, Mesh* mesh, matrix& camera, Light& L) {
	// Combine perspective, camera, and world transformations for the mesh
	matrix p = renderer.perspective * camera * mesh->world;

	unsigned int numCpuCores = std::thread::hardware_concurrency();
	std::vector<std::thread> threads(numCpuCores);

	size_t totalTriangles = mesh->triangles.size();
	size_t work = totalTriangles / numCpuCores;

	for (unsigned int i = 0; i < numCpuCores - 1; i++) {
		threads[i] = std::thread(&renderNew, std::ref(renderer), mesh,
			std::ref(camera), std::ref(L), std::ref(p),
			i * work, (i + 1) * work);
	}

	threads[numCpuCores - 1] = std::thread(&renderNew, std::ref(renderer), mesh,
		std::ref(camera), std::ref(L), std::ref(p),
		(numCpuCores - 1) * work, totalTriangles);

	for (auto& t : threads)
		t.join();
}
*/


// Main rendering function that processes a mesh, transforms its vertices, applies lighting, and draws triangles on the canvas.
// Input Variables:
// - renderer: The Renderer object used for drawing.
// - mesh: Pointer to the Mesh object containing vertices and triangles to render.
// - camera: Matrix representing the camera's transformation.
// - L: Light object representing the lighting parameters.
void render(Renderer& renderer, Mesh* mesh, matrix& camera, Light& L) {
	// Combine perspective, camera, and world transformations for the mesh
	matrix p = renderer.perspective * camera * mesh->world;

	// Iterate through all triangles in the mesh
	for (triIndices& ind : mesh->triangles) {
		Vertex t[3]; // Temporary array to store transformed triangle vertices

		// Transform each vertex of the triangle
		for (unsigned int i = 0; i < 3; i++) {
			t[i].p = p * mesh->vertices[ind.v[i]].p; // Apply transformations
			t[i].p.divideW(); // Perspective division to normalize coordinates

			// Transform normals into world space for accurate lighting
			// no need for perspective correction as no shearing or non-uniform scaling
			t[i].normal = mesh->world * mesh->vertices[ind.v[i]].normal;
			t[i].normal.normalise();

			// Map normalized device coordinates to screen space
			t[i].p[0] = (t[i].p[0] + 1.f) * 0.5f * static_cast<float>(renderer.canvas.getWidth());
			t[i].p[1] = (t[i].p[1] + 1.f) * 0.5f * static_cast<float>(renderer.canvas.getHeight());
			t[i].p[1] = renderer.canvas.getHeight() - t[i].p[1]; // Invert y-axis

			// Copy vertex colours
			t[i].rgb = mesh->vertices[ind.v[i]].rgb;
		}

		// Clip triangles with Z-values outside [-1, 1]
		if (fabs(t[0].p[2]) > 1.0f || fabs(t[1].p[2]) > 1.0f || fabs(t[2].p[2]) > 1.0f) continue;

		// Create a triangle object and render it
		triangle tri(t[0], t[1], t[2]);
		tri.draw(renderer, L, mesh->ka, mesh->kd);
	}
}


// Test scene function to demonstrate rendering with user-controlled transformations
// No input variables
void sceneTest() {
	Renderer renderer;
	// create light source {direction, diffuse intensity, ambient intensity}
	Light L{ vec4(0.f, 1.f, 1.f, 0.f), colour(1.0f, 1.0f, 1.0f), colour(0.1f, 0.1f, 0.1f) };
	// camera is just a matrix
	matrix camera = matrix::makeIdentity(); // Initialize the camera with identity matrix

	bool running = true; // Main loop control variable

	std::vector<Mesh*> scene; // Vector to store scene objects

	// Create a sphere and a rectangle mesh
	Mesh mesh = Mesh::makeSphere(1.0f, 10, 20);
	//Mesh mesh2 = Mesh::makeRectangle(-2, -1, 2, 1);

	// add meshes to scene
	scene.push_back(&mesh);
	// scene.push_back(&mesh2); 

	float x = 0.0f, y = 0.0f, z = -4.0f; // Initial translation parameters
	mesh.world = matrix::makeTranslation(x, y, z);
	//mesh2.world = matrix::makeTranslation(x, y, z) * matrix::makeRotateX(0.01f);

	// Main rendering loop
	while (running) {
		renderer.canvas.checkInput(); // Handle user input
		renderer.clear(); // Clear the canvas for the next frame

		// Apply transformations to the meshes
	 //   mesh2.world = matrix::makeTranslation(x, y, z) * matrix::makeRotateX(0.01f);
		mesh.world = matrix::makeTranslation(x, y, z);

		// Handle user inputs for transformations
		if (renderer.canvas.keyPressed(VK_ESCAPE)) break;
		if (renderer.canvas.keyPressed('A')) x += -0.1f;
		if (renderer.canvas.keyPressed('D')) x += 0.1f;
		if (renderer.canvas.keyPressed('W')) y += 0.1f;
		if (renderer.canvas.keyPressed('S')) y += -0.1f;
		if (renderer.canvas.keyPressed('Q')) z += 0.1f;
		if (renderer.canvas.keyPressed('E')) z += -0.1f;

		// Render each object in the scene
		for (auto& m : scene)
			render(renderer, m, camera, L);

		renderer.present(); // Display the rendered frame
	}
}

// Function to render a scene with multiple objects and dynamic transformations
// No input variables
void scene1() {
	Renderer renderer;
	matrix camera;
	Light L{ vec4(0.f, 1.f, 1.f, 0.f), colour(1.0f, 1.0f, 1.0f), colour(0.1f, 0.1f, 0.1f) };

	bool running = true;

	std::vector<Mesh*> scene;

	// Create a scene of 40 cubes with random rotations
	for (unsigned int i = 0; i < 20; i++) {
		Mesh* m = new Mesh();
		*m = Mesh::makeCube(1.f);
		m->world = matrix::makeTranslation(-2.0f, 0.0f, (-3 * static_cast<float>(i))) * makeRandomRotation();
		scene.push_back(m);
		m = new Mesh();
		*m = Mesh::makeCube(1.f);
		m->world = matrix::makeTranslation(2.0f, 0.0f, (-3 * static_cast<float>(i))) * makeRandomRotation();
		scene.push_back(m);
	}

	float zoffset = 8.0f; // Initial camera Z-offset
	float step = -0.1f;  // Step size for camera movement

	auto start = std::chrono::high_resolution_clock::now();
	std::chrono::time_point<std::chrono::high_resolution_clock> end;
	int cycle = 0;

	// Main rendering loop
	while (running) {
		renderer.canvas.checkInput();
		renderer.clear();

		camera = matrix::makeTranslation(0, 0, -zoffset); // Update camera position

		// Rotate the first two cubes in the scene
		scene[0]->world = scene[0]->world * matrix::makeRotateXYZ(0.1f, 0.1f, 0.0f);
		scene[1]->world = scene[1]->world * matrix::makeRotateXYZ(0.0f, 0.1f, 0.2f);

		if (renderer.canvas.keyPressed(VK_ESCAPE)) break;

		zoffset += step;
		if (zoffset < -60.f || zoffset > 8.f) {
			step *= -1.f;
			if (++cycle % 2 == 0) {
				end = std::chrono::high_resolution_clock::now();
				std::cout << std::chrono::duration<double, std::milli>(end - start).count() << std::endl;
				start = std::chrono::high_resolution_clock::now();
			}
		}


		//for (auto& m : scene)
		//	render(renderer, m, camera, L);

		//unsigned int numCpuCores = std::thread::hardware_concurrency();
		//std::vector<std::thread> threads(numCpuCores);

		//size_t numMesh = scene.size();
		//size_t remain = numMesh % numCpuCores;
		//size_t work = numMesh / numCpuCores;

		//for (unsigned int i = 0; i < numCpuCores - 1; i++) {
		//	threads[i] = std::thread(
		//		[&renderer, &camera, &L, &scene](int start, int end) {
		//			for (unsigned int j = start; j < end; ++j) {
		//				render(renderer, scene[j], camera, L);
		//			}
		//		},
		//		i * work, (i + 1) * work
		//	);
		//}

		//threads[numCpuCores - 1] = std::thread(
		//	[&renderer, &camera, &L, &scene](int start, int end) {
		//		for (unsigned int j = start; j < end; ++j) {
		//			render(renderer, scene[j], camera, L);
		//		}
		//	},
		//	(numCpuCores - 1) * work, numMesh
		//);

		//for (auto& t : threads)
		//	t.join();

		unsigned int numCpuCores = std::thread::hardware_concurrency();
		std::vector<std::thread> threads(numCpuCores);

		size_t numMesh = scene.size();
		size_t remain = numMesh % numCpuCores;
		size_t work = numMesh / numCpuCores;
		size_t startIndex = 0;

		for (unsigned int i = 0; i < numCpuCores; i++) {
			size_t endIndex = startIndex + work + (i < remain ? 1 : 0);
			threads[i] = std::thread(
				[&renderer, &camera, &L, &scene](int start, int end) {
					for (unsigned int j = start; j < end; ++j) {
						render(renderer, scene[j], camera, L);
					}
				},
				startIndex, endIndex
			);
			startIndex = endIndex;
		}

		for (auto& t : threads)
			t.join();

		renderer.present();
	}

	for (auto& m : scene)
		delete m;
}

// Scene with a grid of cubes and a moving sphere
// No input variables
void scene2() {
	Renderer renderer;
	matrix camera = matrix::makeIdentity();
	Light L{ vec4(0.f, 1.f, 1.f, 0.f), colour(1.0f, 1.0f, 1.0f), colour(0.1f, 0.1f, 0.1f) };

	std::vector<Mesh*> scene;

	struct rRot { float x; float y; float z; }; // Structure to store random rotation parameters
	std::vector<rRot> rotations;

	RandomNumberGenerator& rng = RandomNumberGenerator::getInstance();

	// Create a grid of cubes with random rotations
	for (unsigned int y = 0; y < 6; y++) {
		for (unsigned int x = 0; x < 8; x++) {
			Mesh* m = new Mesh();
			*m = Mesh::makeCube(1.f);
			scene.push_back(m);
			m->world = matrix::makeTranslation(-7.0f + (static_cast<float>(x) * 2.f), 5.0f - (static_cast<float>(y) * 2.f), -8.f);
			rRot r{ rng.getRandomFloat(-.1f, .1f), rng.getRandomFloat(-.1f, .1f), rng.getRandomFloat(-.1f, .1f) };
			rotations.push_back(r);
		}
	}

	// Create a sphere and add it to the scene
	Mesh* sphere = new Mesh();
	*sphere = Mesh::makeSphere(1.0f, 10, 20);
	scene.push_back(sphere);
	float sphereOffset = -6.f;
	float sphereStep = 0.1f;
	sphere->world = matrix::makeTranslation(sphereOffset, 0.f, -6.f);

	auto start = std::chrono::high_resolution_clock::now();
	std::chrono::time_point<std::chrono::high_resolution_clock> end;
	int cycle = 0;

	bool running = true;
	while (running) {
		renderer.canvas.checkInput();
		renderer.clear();

		// Rotate each cube in the grid
		for (unsigned int i = 0; i < rotations.size(); i++)
			scene[i]->world = scene[i]->world * matrix::makeRotateXYZ(rotations[i].x, rotations[i].y, rotations[i].z);

		// Move the sphere back and forth
		sphereOffset += sphereStep;
		sphere->world = matrix::makeTranslation(sphereOffset, 0.f, -6.f);
		if (sphereOffset > 6.0f || sphereOffset < -6.0f) {
			sphereStep *= -1.f;
			if (++cycle % 2 == 0) {
				end = std::chrono::high_resolution_clock::now();
				std::cout << std::chrono::duration<double, std::milli>(end - start).count() << std::endl;
				start = std::chrono::high_resolution_clock::now();
			}
		}

		if (renderer.canvas.keyPressed(VK_ESCAPE)) break;

		//for (auto& m : scene)
		//	render(renderer, m, camera, L);

		unsigned int numCpuCores = std::thread::hardware_concurrency();
		std::vector<std::thread> threads(numCpuCores);

		size_t numMesh = scene.size();
		size_t remain = numMesh % numCpuCores;
		size_t work = numMesh / numCpuCores;
		size_t startIndex = 0;

		for (unsigned int i = 0; i < numCpuCores; i++) {
			size_t endIndex = startIndex + work + (i < remain ? 1 : 0);
			threads[i] = std::thread(
				[&renderer, &camera, &L, &scene](int start, int end) {
					for (unsigned int j = start; j < end; ++j) {
						render(renderer, scene[j], camera, L);
					}
				},
				startIndex, endIndex
			);
			startIndex = endIndex;
		}

		for (auto& t : threads)
			t.join();


		renderer.present();
	}

	for (auto& m : scene)
		delete m;
}

void scene3() {
	Renderer renderer;
	// set camera
	matrix camera = matrix::makeTranslation(0, 0, -25);
	// set light
	Light L{ vec4(0.f, 1.f, 1.f, 0.f), colour(1.0f, 1.0f, 1.0f), colour(0.1f, 0.1f, 0.1f) };

	bool running = true;

	std::vector<Mesh*> scene;

	// set central sphere
	Mesh* sphere = new Mesh();
	*sphere = Mesh::makeSphere(5, 12, 24);
	scene.push_back(sphere);
	sphere->world = matrix::makeIdentity();


	// set surronding cubes with layers
	const int cubesPerLayer = 200;
	const int layers = 5;
	const int totalCubes = cubesPerLayer * layers;

	// inner layer has same radius as the sphere so that part would be hide
	// outter layer goes behind the camera and comes back
	float layerR[layers] = { 5,10,15,20,25 };

	std::vector<matrix> cubeTranslationMatrix;
	std::vector<matrix> cubeRotationMatrix;
	std::vector<vec4> cubeRotationSpeed;

	RandomNumberGenerator& rng = RandomNumberGenerator::getInstance();

	for (int l = 0; l < layers; l++) {
		float currentRadius = layerR[l];
		for (int i = 0; i < cubesPerLayer; i++) {
			// evenlly distribute the spheres
			float iFloat = static_cast<float>(i);
			float n = static_cast<float>(cubesPerLayer);
			float phi = acos(1 - 2 * (iFloat + 0.5f) / n);
			float theta = M_PI * (1 + sqrt(5)) * iFloat;
			float x = currentRadius * sin(phi) * cos(theta);
			float y = currentRadius * sin(phi) * sin(theta);
			float z = currentRadius * cos(phi);

			// set transition matrix
			matrix T = matrix::makeTranslation(x, y, z);
			cubeTranslationMatrix.push_back(T);
			cubeRotationMatrix.push_back(matrix::makeIdentity());

			// radom spinning speed
			float rx = rng.getRandomFloat(0.01f, 0.05f) * (rng.getRandomInt(0, 1) ? 1.f : -1.f);
			float ry = rng.getRandomFloat(0.01f, 0.05f) * (rng.getRandomInt(0, 1) ? 1.f : -1.f);
			float rz = rng.getRandomFloat(0.01f, 0.05f) * (rng.getRandomInt(0, 1) ? 1.f : -1.f);
			cubeRotationSpeed.push_back(vec4(rx, ry, rz));

			// create cube meshes
			Mesh* cube = new Mesh();
			*cube = Mesh::makeCube(1.f);
			scene.push_back(cube);

			cube->world = T;
		}
	}

	// general rotate
	float generalRotateAngle = 0.0f;
	float generalRotatationSpeed = 0.01f;
	matrix generalRotationMatrix = matrix::makeRotateY(generalRotateAngle);

	// timing for a complete general rotation
	auto start = std::chrono::high_resolution_clock::now();
	int cycle = 0;

	// main loop
	while (running) {
		renderer.canvas.checkInput();
		renderer.clear();

		// update general angle
		generalRotateAngle += generalRotatationSpeed;
		if (generalRotateAngle >= 2 * M_PI) {
			generalRotateAngle -= 2 * M_PI;
			cycle++;
			auto end = std::chrono::high_resolution_clock::now();
			std::cout << std::chrono::duration<double, std::milli>(end - start).count() << std::endl;
			start = std::chrono::high_resolution_clock::now();
		}
		generalRotationMatrix = matrix::makeRotateY(generalRotateAngle);

		// update sphere spinning
		sphere->world = matrix::makeRotateXYZ(0.02f, 0.03f, 0.04f) * sphere->world;

		// update cube spinning then general

		for (int i = 1; i < totalCubes; i++) {
			cubeRotationMatrix[i] = cubeRotationMatrix[i] *
				matrix::makeRotateXYZ(cubeRotationSpeed[i][0], cubeRotationSpeed[i][1], cubeRotationSpeed[i][2]);

			scene[i]->world = generalRotationMatrix * cubeTranslationMatrix[i] * cubeRotationMatrix[i];
		}


		// render all
		//for (auto& m : scene)
		//	render(renderer, m, camera, L);

		//unsigned int numCpuCores = std::thread::hardware_concurrency();
		unsigned int numCpuCores = 8;

		std::vector<std::thread> threads(numCpuCores);

		size_t numMesh = scene.size();
		size_t remain = numMesh % numCpuCores;
		size_t work = numMesh / numCpuCores;
		size_t startIndex = 0;

		for (unsigned int i = 0; i < numCpuCores; i++) {
			size_t endIndex = startIndex + work + (i < remain ? 1 : 0);
			threads[i] = std::thread(
				[&renderer, &camera, &L, &scene](int start, int end) {
					for (unsigned int j = start; j < end; ++j) {
						render(renderer, scene[j], camera, L);
					}
				},
				startIndex, endIndex
			);
			startIndex = endIndex;
		}

		for (auto& t : threads)
			t.join();


		if (renderer.canvas.keyPressed(VK_ESCAPE)) break;

		renderer.present();
	}

	for (auto& m : scene)
		delete m;

}

// Entry point of the application
// No input variables
int main() {
	// Uncomment the desired scene function to run
	//scene1();
	scene2();
	//scene3();

	//sceneTest(); 


	return 0;
}