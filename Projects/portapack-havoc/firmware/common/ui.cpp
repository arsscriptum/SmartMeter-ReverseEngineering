/*
 * Copyright (C) 2014 Jared Boone, ShareBrained Technology, Inc.
 *
 * This file is part of PortaPack.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#include "ui.hpp"
#include "sine_table.hpp"

#include <algorithm>

namespace ui {

// CGA palette
Color term_colors[16] = {
	Color::black(),
	Color::dark_blue(),
	Color::dark_green(),
	Color::dark_cyan(),
	Color::dark_red(),
	Color::dark_magenta(),
	Color::dark_yellow(),
	Color::light_grey(),
	Color::dark_grey(),
	Color::blue(),
	Color::green(),
	Color::cyan(),
	Color::red(),
	Color::magenta(),
	Color::yellow(),
	Color::white()
};

bool Rect::contains(const Point p) const {
	return (p.x() >= left()) && (p.y() >= top()) &&
	       (p.x() < right()) && (p.y() < bottom());
}

Rect Rect::intersect(const Rect& o) const {
	const auto x1 = std::max(left(), o.left());
	const auto x2 = std::min(right(), o.right());
	const auto y1 = std::max(top(), o.top());
	const auto y2 = std::min(bottom(), o.bottom());
	if( (x2 >= x1) && (y2 > y1) ) {
		return { x1, y1, x2 - x1, y2 - y1 };
	} else {
		return { };
	}
}

// TODO: This violates the principle of least surprise!
// This does a union, but that might not be obvious from "+=" syntax.
Rect& Rect::operator+=(const Rect& p) {
	if( is_empty() ) {
		*this = p;
	}
	if( !p.is_empty() ) {
		const auto x1 = std::min(left(), p.left());
		const auto y1 = std::min(top(), p.top());
		_pos = { x1, y1 };
		const auto x2 = std::max(right(), p.right());
		const auto y2 = std::max(bottom(), p.bottom());
		_size = { x2 - x1, y2 - y1 };
	}
	return *this;
}

Rect& Rect::operator+=(const Point& p) {
	_pos += p;
	return *this;
}

Rect& Rect::operator-=(const Point& p) {
	_pos -= p;
	return *this;
}

Point polar_to_point(float angle, uint32_t distance) {
	return Point(sin_f32(DEG_TO_RAD(angle) + (pi / 2)) * distance, -sin_f32(DEG_TO_RAD(angle)) * distance);
}

} /* namespace ui */
