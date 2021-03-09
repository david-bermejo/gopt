#pragma once

#include <random>
#include <vector>
#include <limits>
#include "vector.hpp"
#include <iostream>
#include <sstream>

#ifdef PLOT_PROGRESS
	#include <matplotlibcpp.h>
	namespace plt = matplotlibcpp;
#endif

#include "../meta/nth_template_type.hpp"

namespace gopt
{
	template <typename V, typename T = typename V::type>
	struct Particle
	{
		V position;
		V velocity;
		V best_position;
		T best_score;
	};

	template <typename T>
	struct PSOptions
	{
		T c1 = 1;
		T c2 = 1;
		T w = 0.9;
		bool print_info = true;
	};

	template <typename F, typename T, unsigned int S>
	Vector_t<T, S> particleswarm(F&& f, const Vector_t<T, S>& lb, const Vector_t<T, S>& ub, const unsigned int n_particles = 2*S, const unsigned int max_iter = 1000, const PSOptions<T>& opts = {})
	{
		return particleswarm(f, lb, ub, S, n_particles, max_iter, opts);
	}

	template <typename F, typename V, typename T = typename V::type>
	V particleswarm(F&& f, const V& lb, const V& ub, const unsigned int size, const unsigned int n_particles, const unsigned int max_iter = 1000, const PSOptions<T>& opts = {})
	{
		static constexpr bool is_dynamic = std::is_same_v<std::remove_cvref_t<V>, Vector<T>>;

		std::default_random_engine gen(time(0));
		std::uniform_real_distribution<> dist(0.0, 1.0);

		const T c1 = opts.c1;
		const T c2 = opts.c2;
		const T w = opts.w;
		const unsigned int max_iter_wo_change = 0.3 * max_iter;

		T g_best_score = std::numeric_limits<T>::max();
		V g_best_position;

		if constexpr (is_dynamic)
			g_best_position = V(size);

		std::vector<T> score_array;
		std::vector<T> indices_array;
		score_array.reserve(max_iter);
		indices_array.reserve(max_iter);

		// Initialize random particles
		std::vector<Particle<V>> particles(n_particles);
		if constexpr (is_dynamic)
		{
			for (auto& p : particles)
			{
				p.position = V(size);
				p.velocity = V(size);
				p.best_position = V(size);
			}
		}

		for (auto& p : particles)
		{
			for (int j = 0; j < size; j++)
				p.position[j] = lb[j] + (ub[j] - lb[j]) * dist(gen);

			if constexpr (is_dynamic) {
				p.velocity.fill(0.0);
			} else {
				p.velocity = 0;
			}
			p.best_position = p.position;
			p.best_score = g_best_score;
		}

		unsigned int iter_wo_change = 0;
		for (int k = 0; k < max_iter; k++)
		{
			const T old_score = g_best_score;

			#pragma omp parallel for
			for (int i = 0; i < n_particles; i++)
			{
				auto& p = particles[i];

				const T score = f(p.position);
				if (score < p.best_score)
				{
					p.best_score = score;
					p.best_position = p.position;

					#pragma omp critical
					if (score < g_best_score)
					{
						g_best_score = score;
						g_best_position = p.position;
					}
				}
			}

			if (g_best_score == old_score)
			{
				iter_wo_change++;
				if (iter_wo_change > max_iter_wo_change)
					break;
			}
			else
				iter_wo_change = 0;

			if (opts.print_info)
				std::cout << " Iteration " << k+1 << ": best score = " << g_best_score << std::endl;

			score_array.push_back(g_best_score);
			indices_array.push_back(k+1);

			// Update particles position
			#pragma omp parallel for
			for (int i = 0; i < n_particles; i++)
			{
				auto& p = particles[i];
				p.velocity = w * p.velocity + (c1 * dist(gen)) * (p.best_position - p.position) + (c2 * dist(gen)) * (g_best_position - p.position);
				p.position += p.velocity;

				for (int j = 0; j < size; j++)
				{
					const T space = ub[j] - lb[j];

					if (p.velocity[j] > space)
						p.velocity[j] = 0;

					else if (p.position[j] > ub[j] || p.position[j] < lb[j])
						p.velocity[j] = -p.velocity[j];

					p.position[j] = std::clamp(p.position[j], lb[j], ub[j]);
				}
			}
		}

#ifdef PLOT_PROGRESS
		plt::figure(1);
		plt::plot(indices_array, score_array, "*");
		plt::title(std::string("Best score: ") + std::to_string(g_best_score));
		plt::xlim(indices_array.front(), indices_array.back());
		plt::grid(true);
#endif

		return g_best_position;
	}
}