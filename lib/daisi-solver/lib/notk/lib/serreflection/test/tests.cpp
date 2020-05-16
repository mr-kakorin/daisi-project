
#include <fstream>
#include <gtest/gtest.h>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <typeinfo>

#include <serreflection/read_json.hpp>
#include <serreflection/write_json.hpp>

#include "test_structures_1.h"
#include "test_structures_2.h"

class TEST_SERREFLECTION : public ::testing::Test
{
  protected:
    virtual void SetUp()
    {
    }
    std::string json_test_1 =
        "{\"v\":\"V1\",\"a\":1,\"b\":[1.1,2.3],\"c\": 1.4,\"e\": true,\"inner\":{\"a\":3, "
        "\"str\":\"str\",\"b\": "
        "[0.1,1.2e17]},\"inners\":[{\"a\":7, "
        "\"b\": "
        "[0.2, "
        "-1.009e15, 3], \"str\":\"str1\"},{\"str\":\"str2\", \"a\":0, \"b\": [0, 5.2]} ]}";

    std::string json_test_2 =
        "{\"v\":\"V1\",\"b\":[6,8],\"a\":1,\"inners\":[{\"a\":7, "
        "\"b\": "
        "[0.2, "
        "-1.009e15, 3], \"str\":\"str1\"},{\"b\": [0, 5.2], \"str\":\"str2\",\"a\":0} ], "
        "\"inner\":{\"str\":\"str\", \"a\":3, \"b\": [0.1, 1.2e17]},\"c\": 1.4,\"e\": true}";

    std::string json_test_3 =
        "{\"v\":\"V1\",\"a\":1,\"b\":[15,15],\"c\": 1.4,\"e\": true,\"inner\":{\"a\":20, "
        "\"str\":\"str\",\"b\": "
        "[0.1,1.2e17]},\"inners\":[{\"b\": "
        "[0.2, "
        "-1.009e15, 3], \"str\":\"str1\"},{\"str\":\"str2\", \"a\":0, \"b\": [0, 5.2]} ]}";

    std::string json_test_4 =
        "{\"v\":\"V1\",\"b\":[1.1,2.3],\"a\":1,\"inners\":[{\"a\":7, "
        "\"b\": "
        "[0.2, "
        "-1.009e15, 3], \"str\":\"str1\"},{\"b\": [0, 5.2], \"str\":\"str2\",\"a\":0} ], "
        "\"inner\":{\"str\":\"str\", \"a\":3, \"b\": [0.1,1.2e17]}}";

    std::string json_test_5 =
        "{\"v\":\"V4\",\"b\":[1.1,2.3],\"a\":1,\"inners\":[{\"a\":7, "
        "\"b\": "
        "[0.2, "
        "-1.009e15, 3], \"str\":\"str1\"},{\"b\": [0, 5.2], \"str\":\"str2\",\"a\":0} ], "
        "\"inner\":{\"str\":\"str\", \"a\":3, \"b\": [0.1,1.2e17]}}";

    std::string json_test_6 =
        "{\"b\":[1.1,2.3],\"a\":1,\"inners\":[{\"a\":7, "
        "\"b\": "
        "[0.2, "
        "-1.009e15, 3], \"str\":\"str1\"},{\"b\": [0, 5.2], \"str\":\"str2\",\"a\":0} ], "
        "\"inner\":{\"str\":\"str\", \"a\":3, \"b\": [0.1,1.2e17]}}";
};

template <class T>
void check_inner(const T& actual, const int a_exp, std::vector<double> b_exp,
                 const std::string& str_exp)
{
    EXPECT_EQ(a_exp, actual.a);
    ASSERT_EQ(b_exp.size(), actual.b.size());
    for (size_t i = 0; i < b_exp.size(); i++)
    {
        EXPECT_DOUBLE_EQ(b_exp[i], actual.b[i]);
    }
    EXPECT_STREQ(str_exp.c_str(), actual.str.c_str());
}

template <class T>
void test(const T& actual, double b0_exp, double b1_exp)
{
    EXPECT_EQ(1, actual.a);
    ASSERT_EQ(2, actual.b.size());
    EXPECT_DOUBLE_EQ(b0_exp, actual.b[0]);
    EXPECT_DOUBLE_EQ(b1_exp, actual.b[1]);

    EXPECT_DOUBLE_EQ(1.4, actual.c);
    EXPECT_TRUE(actual.e);

    check_inner(actual.inner, 20, {0.1, 1.2e17}, "str");
    EXPECT_EQ(2, actual.inners.size());
    check_inner(actual.inners[0], 20, {0.2, -1.009e15, 3}, "str1");
    check_inner(actual.inners[1], 20, {0, 5.2}, "str2");
}

TEST_F(TEST_SERREFLECTION, test1)
{
    ASSERT_FALSE(srfl::read_json_string<TestJSON_1>(json_test_1));

    auto res_2 = srfl::read_json_string<TestJSON_1>(json_test_2);

    ASSERT_TRUE(res_2);

    test(*res_2, 6, 8);

    auto str = srfl::write_json_str(*res_2);

    auto res_22 = srfl::read_json_string<TestJSON_1>(str);

    ASSERT_TRUE(res_22);

    test(*res_22, 6, 8);
}
TEST_F(TEST_SERREFLECTION, test2)
{
    ASSERT_FALSE(srfl::read_json_string<TestJSON_1>(json_test_3));
    ASSERT_FALSE(srfl::read_json_string<TestJSON_1>(json_test_4));
}
TEST_F(TEST_SERREFLECTION, test3)
{
    auto res_1 = srfl::read_json_string<test_srfl::TestJSON_2>(json_test_1);

    ASSERT_TRUE(res_1);

    test(*res_1, 4, 4);

    auto str_res_1 = srfl::write_json_str(*res_1);

    auto res_1_j = srfl::read_json_string<test_srfl::TestJSON_2>(str_res_1);

    ASSERT_TRUE(res_1_j);

    test(*res_1_j, 4, 4);

    auto res_2 = srfl::read_json_string<test_srfl::TestJSON_2>(json_test_2);

    ASSERT_TRUE(res_2);

    test(*res_2, 6, 8);

    auto str_res_2 = srfl::write_json_str(*res_2);

    auto res_2_j = srfl::read_json_string<test_srfl::TestJSON_2>(str_res_2);

    ASSERT_TRUE(res_2_j);

    test(*res_2_j, 6, 8);
}

TEST_F(TEST_SERREFLECTION, test4)
{
    auto res_1 = srfl::read_json_string<test_srfl::TestJSON_2>(json_test_3);

    ASSERT_TRUE(res_1);

    EXPECT_EQ(res_1->v, test_srfl::En::V1);

    test(*res_1, 15, 15);

    auto str_res_1 = srfl::write_json_str(*res_1);

    auto res_1_j = srfl::read_json_string<test_srfl::TestJSON_2>(str_res_1);

    ASSERT_TRUE(res_1_j);

    test(*res_1_j, 15, 15);

    auto res_2 = srfl::read_json_string<test_srfl::TestJSON_2>(json_test_4);

    ASSERT_TRUE(res_2);

    test(*res_2, 4, 4);

    auto str_res_2 = srfl::write_json_str(*res_2);

    auto res_2_j = srfl::read_json_string<test_srfl::TestJSON_2>(str_res_2);

    ASSERT_TRUE(res_2_j);

    test(*res_2_j, 4, 4);

    auto res_3 = srfl::read_json_string<test_srfl::TestJSON_2>(json_test_5);

    ASSERT_TRUE(res_3);

    EXPECT_EQ(unsigned(res_3->v), 1);

    auto res_4 = srfl::read_json_string<test_srfl::TestJSON_2>(json_test_6);

    ASSERT_TRUE(res_4);

    EXPECT_EQ(unsigned(res_4->v), 1);
}
