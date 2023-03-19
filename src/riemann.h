/*
 * Copyright (c) 2023 Bryance Oyang
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef RIEMANN_H
#define RIEMANN_H

#include "array.h"
#include "grid.h"

namespace riemann {

void HLLC(const Array<double> &Lprim, const Array<double> &Lcons, const Array<double> &Lw_array,
	const Array<double> &Rprim, const Array<double> &Rcons, const Array<double> &Rw_array,
	Array<double> &J, int dir, int il, int iuf, int jl, int ju);

/*void HLLE(const Array<double> &Lcons, const Array<double> &LJ_array, const Array<double> &Lw_array,
	const Array<double> &Rcons, const Array<double> &RJ_array, const Array<double> &Rw_array,
	Array<double> &J, int dir);*/

} // namespace riemann

#endif /* RIEMANN_H */
