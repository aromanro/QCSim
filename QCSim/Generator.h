#pragma once

#include <vector>
#include <algorithm>

#include "PauliStringXZ.h"

namespace QC {
	namespace Clifford {

		class Generator : public PauliStringXZWithSign {
		public:
			Generator() : PauliStringXZWithSign() {}

			Generator(size_t nQubits) : PauliStringXZWithSign(nQubits) {}

			Generator(const Generator& other) : PauliStringXZWithSign(other) {}

			Generator(Generator&& other) noexcept : PauliStringXZWithSign(std::move(other)) {}

			Generator& operator=(const Generator& other)
			{
				if (this != &other)
				{
					PauliStringXZWithSign::operator=(other);
				}

				return *this;
			}

			Generator& operator=(Generator&& other) noexcept
			{
				if (this != &other)
				{
					PauliStringXZWithSign::operator=(std::move(other));
				}

				return *this;
			}
		};

	}
}
