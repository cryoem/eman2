/*
 *
 * Copyright (c) 2020- Baylor College of Medicine
 * 
 * This software is issued under a joint BSD/GNU license. You may use the
 * source code in this file under either license. However, note that the
 * complete EMAN2 and SPARX software packages have some GPL dependencies,
 * so you are responsible for compliance with the licenses of these packages
 * if you opt to use BSD licensing. The warranty disclaimer below holds
 * in either instance.
 * 
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 * 
 * */

#include "io/eerio.h"

using namespace EMAN;


#undef NDEBUG
#include <cassert>


#include "io/eerio.h"

void test_bit_stream() {
	uint64_t a1 = ~0;

	BitStream<uint64_t> is1(&a1);
	assert(is1.get_bits(3) == 7);
	assert(is1.get_bits(9) == (1 << 9) - 1);
	assert(is1.get_bits(1) == 1);

	uint64_t a2 = 1<<1 | 1<<3 | 1<<5 | 1<<7;
	BitStream<uint64_t> is2(&a2);

	assert(is2.get_bits(2) == 2);
	assert(is2.get_bits(3) == 2);
	assert(is2.get_bits(3) == 5);

	uint8_t a3 = 0b10101010;
	uint8_t b3 = 0b10011001;
	uint8_t ab3[] = {a3, b3};

	BitStream<uint64_t> is3(reinterpret_cast<uint64_t *>(ab3));
	assert(is3.get_bits(2) == 2);
	assert(is3.get_bits(6) == 42);
	assert(is3.get_bits(3) == 1);

	uint64_t ab4[] = {a3, b3};
	BitStream<uint64_t> is4(ab4);

	assert(is4.get_bits(20) == a3);
	assert(is4.get_bits(43) == 0);
	assert(is4.get_bits(20) == 0b100110010);
}

uint8_t a5 = 0b10101010;
uint8_t b5 = 0b11011001;
uint8_t c5 = 0b10011111;
uint8_t d5 = 0b10011001;
uint8_t ab5[] = {a5, b5, c5, d5};

void test_bit_reader() {
	BitStream<uint8_t> is5(ab5);
	Rle<7, uint8_t> rle;

	is5 >> rle;
	assert(rle == 42);

	is5 >> rle;
	assert(rle == 0b110011);

	is5 >> rle;
	assert(rle == 0b1001100 + (1<<7) - 1);
}

void test_eer_sub_pix() {
	BitStream<uint8_t> is6(ab5);
	SubPix<4, uint8_t> sub_pix;
	Rle<7, uint8_t> rle6;

	is6 >> sub_pix;
	assert(sub_pix == 10);

	is6 >> rle6 >> sub_pix;
	assert(sub_pix == 11);
}

void test_eer_rle_no_overflow() {
	typedef uint8_t BuffWord;
	auto num_bits = sizeof(BuffWord) * 8;
	
	BuffWord a = 0b11111111;
	
	BitStream<BuffWord> is1(&a);
	BitReader<1, false, BuffWord> rle1;
	
	for(int i=0; i<num_bits; i++) {
		is1 >> rle1;
		assert(rle1 == 1);
	}

	BitStream<BuffWord> is2(&a);
	BitReader<2, false, BuffWord> rle2;

	for(int i=0; i<num_bits/2; i++) {
		is2 >> rle2;
		assert(rle2 == 0b11);
	}

	BitStream<BuffWord> is4(&a);
	BitReader<4, false, BuffWord> rle4;

	for(int i=0; i<num_bits/4; i++) {
		is4 >> rle4;
		assert(rle4 == 0b1111);
	}

	BitStream<BuffWord> is3(&a);
	BitReader<3, false, BuffWord> rle3;

	is3>>rle3;
	assert(rle3 == 0b111);
	is3>>rle3;
	assert(rle3 == 0b111);
	is3>>rle3;
	assert(rle3 == 0b11);
}

int main()
{
	test_bit_stream();
	test_bit_reader();
	test_eer_sub_pix();
	test_eer_rle_no_overflow();

	return 0;
}
