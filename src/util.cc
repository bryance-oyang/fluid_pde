/*
 * Copyright (c) 2023 Bryance Oyang
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include <cmath>
#include "macro.hh"

namespace util {

number fmin3(number a, number b, number c)
{
	return fmin(fmin(a,b),c);
}

number fmin4(number a, number b, number c, number d)
{
	return fmin(fmin(fmin(a,b),c),d);
}

}
