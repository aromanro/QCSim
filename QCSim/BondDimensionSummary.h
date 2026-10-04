#pragma once

#include <algorithm>
#include <cstddef>
#include <vector>

namespace QC
{
namespace TensorNetworks
{
// Optional maximum tree. Local dimension changes allocate nothing; unchanged
// dimensions cost one comparison. Disabled simulators skip tree accesses.
template <typename IndexType> class BondDimensionSummary
{
  public:
    template <class Lambdas> void Enable(bool enabled, const Lambdas &lambdas)
    {
        if (!enabled)
        {
            tree.clear();
            leafCount = 0;
            return;
        }
        leafCount = 1;
        while (leafCount < lambdas.size())
            leafCount *= 2;
        tree.assign(2 * leafCount, IndexType(1));
        Refresh(lambdas);
    }

    template <class Lambdas> void Refresh(const Lambdas &lambdas)
    {
        if (tree.empty())
            return;
        if (lambdas.size() > leafCount)
        {
            Enable(true, lambdas);
            return;
        }
        std::fill(tree.begin() + leafCount, tree.end(), IndexType(1));
        for (size_t i = 0; i < lambdas.size(); ++i)
            tree[leafCount + i] = lambdas[i].size();
        for (size_t i = leafCount - 1; i != 0; --i)
            tree[i] = std::max(tree[2 * i], tree[2 * i + 1]);
    }

    void Update(size_t bond, IndexType dimension) noexcept
    {
        if (tree.empty())
            return;
        size_t i = leafCount + bond;
        if (tree[i] == dimension)
            return;
        tree[i] = dimension;
        while ((i /= 2) != 0)
        {
            const IndexType value = std::max(tree[2 * i], tree[2 * i + 1]);
            if (tree[i] == value)
                break;
            tree[i] = value;
        }
    }

    template <class Lambdas> IndexType Maximum(const Lambdas &lambdas) const
    {
        if (!tree.empty())
            return tree[1];
        IndexType maximum = 1;
        for (const auto &lambda : lambdas)
            maximum = std::max(maximum, lambda.size());
        return maximum;
    }

  private:
    size_t leafCount = 0;
    std::vector<IndexType> tree;
};
} // namespace TensorNetworks
} // namespace QC
