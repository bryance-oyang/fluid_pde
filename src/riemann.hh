/*
 * Copyright (c) 2023 Bryance Oyang
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef RIEMANN_H
#define RIEMANN_H

#include "array.hh"
#include "grid.hh"

namespace riemann {

void HLLC(const Array<number> &Lprim, const Array<number> &Lcons, const Array<number> &Lw_array,
	const Array<number> &Rprim, const Array<number> &Rcons, const Array<number> &Rw_array,
	Array<number> &J, int dir, int il, int iuf, int jl, int ju);

void HLLE(const Array<number> &Lcons, const Array<number> &LJ_array, const Array<number> &Lw_array,
	const Array<number> &Rcons, const Array<number> &RJ_array, const Array<number> &Rw_array,
	Array<number> &J, int dir, int il, int iuf, int jl, int ju);

} // namespace riemann

#endif /* RIEMANN_H */
