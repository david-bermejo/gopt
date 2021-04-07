#pragma once

#include <vector>
#include <tuple>
#include <random>
#include <numeric>
#include "vector.hpp"

namespace gopt
{

    template <typename T>
    std::vector<unsigned int> tournament_selection(const std::vector<T>& fitness)
    {
        const unsigned int size = fitness.size();

        std::vector<unsigned int> mating_pool(size);
        std::vector<unsigned int> index(size);
        for (int i = 0; i < size; i++)
            index[i] = i;
        
        std::default_random_engine gen(time(0));
        std::shuffle(index.begin(), index.end(), gen);

        for (int i = 0; i < size-1; i++)
        {
            const unsigned int id1 = index[i];
            const unsigned int id2 = index[i+1];
            mating_pool[i] = (fitness[id1] < fitness[id2]) ? id1 : id2;
        }

        const unsigned int id1 = index[size-1];
        const unsigned int id2 = index[0];
        mating_pool[size-1] = (fitness[id1] < fitness[id2]) ? id1 : id2;

        return mating_pool;
    }

    template <typename V, typename T = typename V::type>
    std::vector<V> crossover_SBX(const std::vector<V>& parents, const T p_cross, const T eta_cross, const V& lb, const V& ub)
    {
        constexpr bool is_dynamic = std::is_same_v<std::remove_cvref_t<V>, Vector<T>>;

        const unsigned int Np = parents.size();
        const unsigned int S = parents[0].size();

        std::vector<V> offsprings(Np);
        std::vector<unsigned int> indices(Np);
        for (int i = 0; i < Np; i++)
        {
            indices[i] = i;
            if constexpr (is_dynamic)
                offsprings[i] = V(S);
        }
        
        std::default_random_engine gen(time(0));
        std::uniform_real_distribution<T> dist(0, 1);
        std::shuffle(indices.begin(), indices.end(), gen);

        for (int i = 0; i < Np; i += 2)
        {
            const unsigned int id1 = indices[i];
            const unsigned int id2 = indices[i+1];

            T u = dist(gen);

            if (u < p_cross)
            {
                for (int j = 0; j < S; j++)
                {
                    u = dist(gen);

                    T beta;
                    if (u <= 0.5)
                        beta = std::pow(2 * u, 1.0 / (eta_cross + 1));
                    else
                        beta = std::pow(0.5/(1 - u), 1.0 / (eta_cross + 1));
                    
                    offsprings[i][j] = 0.5 * ((1 + beta) * parents[id1][j] + (1 - beta) * parents[id2][j]);
                    offsprings[i+1][j] = 0.5 * ((1 - beta) * parents[id1][j] + (1 + beta) * parents[id2][j]);
                }

                for (int j = 0; j < S; j++)
                {
                    offsprings[i][j] = std::clamp(offsprings[i][j], lb[j], ub[j]);
                    offsprings[i+1][j] = std::clamp(offsprings[i+1][j], lb[j], ub[j]);
                }
            }
            else
            {
                offsprings[i] = parents[id1];
                offsprings[i+1] = parents[id2];
            }
        }

        return offsprings;
    }

    template <typename V, typename T = typename V::type>
    std::vector<V> mutation_poly(std::vector<V> offsprings, const T p_mut, const T eta_mut, const V& lb, const V& ub)
    {
        const unsigned int Np = offsprings.size();
        const unsigned int S = offsprings[0].size();

        std::default_random_engine gen(time(0));
        std::uniform_real_distribution<double> dist(0, 1);

        for (int i = 0; i < Np; i++)
        {
            T u = dist(gen);

            if (u < p_mut)
            {
                for (int j = 0; j < S; j++)
                {
                    u = dist(gen);

                    T delta;
                    if (u < 0.5)
                        delta = std::pow(2*u, 1.0 / (eta_mut + 1)) - 1;
                    else
                        delta = 1 - std::pow(2*(1-u), 1.0 / (eta_mut + 1));
                    
                    offsprings[i][j] += (ub[j] - lb[j]) * delta;
                    offsprings[i][j] = std::clamp(offsprings[i][j], lb[j], ub[j]);
                }
            }
        }

        return offsprings;
    }

    template <typename T>
    std::vector<size_t> sort_indices(const std::vector<T>& v)
    {
        std::vector<size_t> idx(v.size());
        std::iota(idx.begin(), idx.end(), 0);

        std::stable_sort(idx.begin(), idx.end(),
            [&](size_t i1, size_t i2) { return v[i1] < v[i2]; });
        
        return idx;
    }

    template <typename F, typename V, typename T = typename V::type>
    std::tuple<V, T> ga(F&& f, const std::vector<V>& x0, const V& lb, const V& ub, const T p_mut, const T p_cross, const T eta_mut, const T eta_cross, const unsigned int max_iter = 1000)
    {
        // Population size
        const unsigned int Np = x0.size();
        const unsigned int S = x0[0].size();

        // Initial population
        std::vector<V> population = x0;
        std::vector<T> population_fitness(Np);
        std::vector<T> offspring_fitness(Np);

        // Evaluate fitness function for initial population
        #pragma omp parallel for
        for (int i = 0; i < Np; i++)
            population_fitness[i] = f(population[i]);
        
        std::vector<double> xaxis;
        std::vector<double> yaxis;

        // Ieration loop
        for (int k = 0; k < max_iter; k++)
        {
            // Tournament selection
            const auto indices = tournament_selection(population_fitness);
            std::vector<V> parents(Np);
            for (int i = 0; i < Np; i++)
                parents[i] = population[indices[i]];
            
            // Crossover
            auto offsprings = crossover_SBX(parents, p_cross, eta_cross, lb, ub);

            // Mutation
            offsprings = mutation_poly(offsprings, p_mut, eta_mut, lb, ub);

            #pragma omp parallel for
            for (int j = 0; j < Np; j++)
                offspring_fitness[j] = f(offsprings[j]);
            
            // Combine
            std::vector<V> comb_pop;
            comb_pop.reserve(2 * Np);
            comb_pop.insert(comb_pop.end(), population.begin(), population.end());
            comb_pop.insert(comb_pop.end(), offsprings.begin(), offsprings.end());
            
            std::vector<T> comb_fit;
            comb_fit.reserve(2 * Np);
            comb_fit.insert(comb_fit.end(), population_fitness.begin(), population_fitness.end());
            comb_fit.insert(comb_fit.end(), offspring_fitness.begin(), offspring_fitness.end());

            std::vector<size_t> sorted_indices = sort_indices(comb_fit);

            for (int i = 0; i < Np; i++)
            {
                population[i] = comb_pop[sorted_indices[i]];
                population_fitness[i] = comb_fit[sorted_indices[i]];
            }

            if (true)
            {
                xaxis.push_back(k);
                yaxis.push_back(population_fitness[0]);

                /*plt::plot(xaxis, yaxis, "red");
                plt::xlabel("Iteration");
                plt::ylabel("Cost");
                plt::title(std::string("Score: ") + std::to_string(population_fitness[0]));
                plt::draw();
                plt::pause(0.001);*/
            }
        }

        const size_t index = std::min_element(population_fitness.begin(), population_fitness.end()) - population_fitness.begin();
        return { population[index], population_fitness[index] };
    }
}