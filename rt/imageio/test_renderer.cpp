#include "io/renderer.h"

using namespace EMAN;

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>

struct Tester : public Renderer
{
	using Renderer::getRenderedDataAndRendertrunc;
	using Renderer::rendered_dt;

	using Renderer::renderbits;
	using Renderer::rendermax;
	using Renderer::rendermin;
};

Tester a;

TEMPLATE_TEST_CASE("scale to bits=2, min=0, max=3", "", unsigned char, unsigned short, unsigned int) {
	vector<float> data {-100000, 1, 2, 100000};
	a.renderbits = 2;  // 2^2 = 4
	a.rendermax = 4 - 1; // 3
	a.rendermin = 0;
	auto [vi, count] = a.getRenderedDataAndRendertrunc<TestType>(data.data(), data.size());
	REQUIRE(count == 2);
	REQUIRE(vi == vector<TestType>{0, 1, 2, 3});
}

TEMPLATE_TEST_CASE("scale to bits=3, min=0, max=7", "", unsigned char, unsigned short, unsigned int) {
	vector<float> data{-100000, 1, 2, 100000};
	a.renderbits = 3;  // 2^3 = 8
	a.rendermax = 8 - 1; // 7
	a.rendermin = 0;
	auto [vi, count] = a.getRenderedDataAndRendertrunc<TestType>(data.data(), data.size());
	REQUIRE(count == 2);
	REQUIRE(vi == vector<TestType>{0, 1, 2, 7});
}

TEMPLATE_TEST_CASE("scale to bits=3, min=3, max=10", "", unsigned char, unsigned short, unsigned int) {
	vector<float> data{-100000, 1, 2, 100000};
	a.renderbits = 3;  // 2^3 = 8
	a.rendermax = 10; // 7 + 3
	a.rendermin = 3;  // 0 + 3
	auto [vi1, count1] = a.getRenderedDataAndRendertrunc<TestType>(data.data(), data.size());
	REQUIRE(count1 == 4);
	REQUIRE(vi1 == vector<TestType>{0, 0, 0, 7});
}

TEMPLATE_TEST_CASE("scale to bits=4, min=0, max=15", "", unsigned char, unsigned short, unsigned int) {
	vector<float> data {-100000, 1, 2, 10, 11, 100000};
	a.renderbits = 4;  // 2^4 = 16
	a.rendermax = 16 - 1; // 15
	a.rendermin = 0;
	auto [vi, count] = a.getRenderedDataAndRendertrunc<TestType>(data.data(), data.size());
	REQUIRE(count == 2);
	REQUIRE(vi == vector<TestType>{0, 1, 2, 10, 11, 15});
}

TEMPLATE_TEST_CASE("scale to bits=7, min=0, max=127", "", unsigned char, unsigned short, unsigned int) {
	vector<float> data {-100000, 1, 2, 10, 11, 100000};
	a.renderbits = 7; // 2^7 = 128
	a.rendermax = 128 - 1; // 127
	a.rendermin = 0;
	auto [vi, count] = a.getRenderedDataAndRendertrunc<TestType>(data.data(), data.size());
	REQUIRE(count == 2);
	REQUIRE(vi == vector<TestType>{0, 1, 2, 10, 11, 127});
}

TEMPLATE_TEST_CASE("scale to bits=7, min=42, max=169", "", unsigned char, unsigned short, unsigned int) {
	float sh = 42;
	vector<float> data {-100000, 1+sh, 2+sh, 10+sh, 11+sh, 100000};
	a.renderbits = 7; // 2^7 = 128
	a.rendermax = 128 - 1 + sh;
	a.rendermin = 0       + sh;
	auto [vi, count] = a.getRenderedDataAndRendertrunc<TestType>(data.data(), data.size());
	REQUIRE(count == 2);
	REQUIRE(vi == vector<TestType>{0, 1, 2, 10, 11, 127});
}

TEMPLATE_TEST_CASE("scale to bits=4, min=0, max=15", "", char, short, int) {
	vector<float> data {-100000, 0, 2, 10, 14, 100000};
	a.renderbits = 4;  // 2^4 = 16
	a.rendermax = 16 - 1; // 15
	a.rendermin = 0;
	auto [vi, count] = a.getRenderedDataAndRendertrunc<TestType>(data.data(), data.size());
	REQUIRE(count == 3);
	REQUIRE(vi == vector<TestType>{-8, -8, -6, 2, 6, 7});
}

TEMPLATE_TEST_CASE("scale to bits=6, min=0, max=15", "", char, short, int) {
	vector<float> data {-100000, 0, 2, 10, (73+10)/2 + 1, 100000};
	a.renderbits = 6;  // 2^6 = 64
	a.rendermax = 73;  // 63 + 10
	a.rendermin = 10;  // 0 + 10
	auto [vi, count] = a.getRenderedDataAndRendertrunc<TestType>(data.data(), data.size());
	REQUIRE(count == 5);
	REQUIRE(vi == vector<TestType>{-32, -32, -32, -32, 0, 31});
}

TEST_CASE("return datatype when compression requested") {
	a.renderbits = 16;
	REQUIRE(a.rendered_dt(EMUtil::EM_COMPRESSED, {EMUtil::EM_UCHAR, EMUtil::EM_USHORT, EMUtil::EM_FLOAT}) == EMUtil::EM_USHORT);
	a.renderbits = 8;
	REQUIRE(a.rendered_dt(EMUtil::EM_COMPRESSED, {EMUtil::EM_UCHAR, EMUtil::EM_USHORT, EMUtil::EM_FLOAT}) == EMUtil::EM_UCHAR);
	a.renderbits = 0;
	REQUIRE(a.rendered_dt(EMUtil::EM_COMPRESSED, {EMUtil::EM_UCHAR, EMUtil::EM_USHORT, EMUtil::EM_FLOAT}) == EMUtil::EM_FLOAT);

	a.renderbits = 12;
	REQUIRE(a.rendered_dt(EMUtil::EM_COMPRESSED, {EMUtil::EM_UCHAR, EMUtil::EM_USHORT, EMUtil::EM_FLOAT}) == EMUtil::EM_USHORT);
	a.renderbits = 4;
	REQUIRE(a.rendered_dt(EMUtil::EM_COMPRESSED, {EMUtil::EM_UCHAR, EMUtil::EM_USHORT, EMUtil::EM_FLOAT}) == EMUtil::EM_UCHAR);

	a.renderbits = 16;
	REQUIRE(a.rendered_dt(EMUtil::EM_COMPRESSED, {EMUtil::EM_UCHAR, EMUtil::EM_USHORT, EMUtil::EM_FLOAT}) == EMUtil::EM_USHORT);
}

TEST_CASE("various initializer lists when compression requested") {
	REQUIRE(a.rendered_dt(EMUtil::EM_COMPRESSED, {EMUtil::EM_UCHAR, EMUtil::EM_SHORT}) == EMUtil::EM_SHORT);
	a.renderbits = 8;
	REQUIRE(a.rendered_dt(EMUtil::EM_COMPRESSED, {EMUtil::EM_UCHAR, EMUtil::EM_USHORT}) == EMUtil::EM_UCHAR);
	REQUIRE(a.rendered_dt(EMUtil::EM_COMPRESSED, {EMUtil::EM_CHAR, EMUtil::EM_USHORT}) == EMUtil::EM_CHAR);
	REQUIRE(a.renderbits == 8);

	a.renderbits = 9;
	REQUIRE(a.rendered_dt(EMUtil::EM_COMPRESSED, {EMUtil::EM_CHAR, EMUtil::EM_FLOAT, EMUtil::EM_SHORT}) == EMUtil::EM_SHORT);
	a.renderbits = 4;
	REQUIRE(a.rendered_dt(EMUtil::EM_COMPRESSED, {EMUtil::EM_CHAR, EMUtil::EM_FLOAT, EMUtil::EM_SHORT}) == EMUtil::EM_CHAR);
	REQUIRE(a.rendered_dt(EMUtil::EM_COMPRESSED, {EMUtil::EM_FLOAT, EMUtil::EM_SHORT}) == EMUtil::EM_SHORT);
	a.renderbits = 0;
	REQUIRE(a.rendered_dt(EMUtil::EM_COMPRESSED, {EMUtil::EM_FLOAT, EMUtil::EM_SHORT}) == EMUtil::EM_FLOAT);
}

TEST_CASE("return datatype when no compression requested") {
	REQUIRE(a.rendered_dt(EMUtil::EM_USHORT, {EMUtil::EM_UCHAR, EMUtil::EM_USHORT, EMUtil::EM_FLOAT}) == EMUtil::EM_USHORT);
	REQUIRE(a.rendered_dt(EMUtil::EM_UCHAR, {EMUtil::EM_UCHAR, EMUtil::EM_USHORT, EMUtil::EM_FLOAT}) == EMUtil::EM_UCHAR);
	REQUIRE(a.rendered_dt(EMUtil::EM_FLOAT, {EMUtil::EM_UCHAR, EMUtil::EM_USHORT, EMUtil::EM_FLOAT}) == EMUtil::EM_FLOAT);
	REQUIRE(a.rendered_dt(EMUtil::EM_INT, {EMUtil::EM_UCHAR, EMUtil::EM_USHORT, EMUtil::EM_FLOAT}) == EMUtil::EM_INT);
	REQUIRE(a.rendered_dt(EMUtil::EM_UINT, {EMUtil::EM_UCHAR, EMUtil::EM_USHORT, EMUtil::EM_FLOAT}) == EMUtil::EM_UINT);
}

TEST_CASE("prefer unsigned type over signed") {
	a.renderbits = 8;
	REQUIRE(a.rendered_dt(EMUtil::EM_COMPRESSED, {EMUtil::EM_CHAR, EMUtil::EM_UCHAR}) == EMUtil::EM_UCHAR);
	REQUIRE(a.rendered_dt(EMUtil::EM_COMPRESSED, {EMUtil::EM_UCHAR, EMUtil::EM_CHAR}) == EMUtil::EM_UCHAR);
	REQUIRE(a.renderbits == 8);
	REQUIRE(a.rendered_dt(EMUtil::EM_COMPRESSED, {EMUtil::EM_SHORT, EMUtil::EM_USHORT}) == EMUtil::EM_USHORT);
	REQUIRE(a.rendered_dt(EMUtil::EM_COMPRESSED, {EMUtil::EM_USHORT, EMUtil::EM_SHORT}) == EMUtil::EM_USHORT);
	REQUIRE(a.renderbits == 8);

	a.renderbits = 12;
	REQUIRE(a.rendered_dt(EMUtil::EM_COMPRESSED, {EMUtil::EM_SHORT, EMUtil::EM_USHORT}) == EMUtil::EM_USHORT);
	REQUIRE(a.rendered_dt(EMUtil::EM_COMPRESSED, {EMUtil::EM_USHORT, EMUtil::EM_SHORT}) == EMUtil::EM_USHORT);
	REQUIRE(a.renderbits == 12);

	a.renderbits = 16;
	REQUIRE(a.rendered_dt(EMUtil::EM_COMPRESSED, {EMUtil::EM_SHORT, EMUtil::EM_USHORT}) == EMUtil::EM_USHORT);
	REQUIRE(a.renderbits == 16);
	REQUIRE(a.rendered_dt(EMUtil::EM_COMPRESSED, {EMUtil::EM_USHORT, EMUtil::EM_SHORT}) == EMUtil::EM_USHORT);
	REQUIRE(a.renderbits == 16);

	a.renderbits = 4;
	REQUIRE(a.rendered_dt(EMUtil::EM_COMPRESSED, {EMUtil::EM_CHAR, EMUtil::EM_UCHAR}) == EMUtil::EM_UCHAR);
	REQUIRE(a.renderbits == 4);
	REQUIRE(a.rendered_dt(EMUtil::EM_COMPRESSED, {EMUtil::EM_UCHAR, EMUtil::EM_CHAR}) == EMUtil::EM_UCHAR);
	REQUIRE(a.renderbits == 4);
}
