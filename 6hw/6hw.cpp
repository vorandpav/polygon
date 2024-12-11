#include <bits/stdc++.h>

using namespace std;

class Solution {
public:
    bool check(const string &s, int k, const string word) {
        // проверка, что в строке s есть k подпоследовательностей word
        int word_index = 0;
        for (char c: s) {
            if (c == word[word_index % word.size()])
                word_index++;
            if (word_index == k * word.size())
                return true;
        }
        return false;

        // O(n) - время
    }

    string longestSubsequenceRepeatedK(string s, int k) {
        // n = s.size()
        vector<int> letter(26);
        // массив количеств букв в строке s

        for (char c: s)
            letter[c - 'a']++;
        // O(n)

        for (int &amount: letter)
            amount /= k;
        // считаем сколько букв можно использовать, чтобы они вместе могли встретиться k раз

        set<pair<string, vector<int> > > words;
        // множество пар: {слово, массив неиспользованных букв}

        words.insert({"", letter});

        string result = "";

        while (!words.empty()) {
            set<pair<string, vector<int> > > new_words;
            // новые более длинные слова

            for (auto word = words.rbegin(); word != words.rend(); ++word) {
                // в порядке лексикографического убывания

                for (int i = 25; i > -1; --i) {
                    if (word->second[i] > 0) {
                        if (check(s, k, word->first + (char) ('a' + i))) {
                            // если новое слово встретилось k раз в строке s

                            vector<int> new_letter = word->second;
                            --new_letter[i];
                            new_words.insert({word->first + (char) ('a' + i), new_letter});

                            if (word->first.size() + 1 > result.size())
                                result = word->first + (char) ('a' + i);
                        }
                    }
                }
            }
            // O(n) - время
            words = new_words;
        }
        // перебор длин может дойти до n/k
        // O(n * n/k) - время

        return result;
    }

    // O(n * n/k) - время
};

int main() {
    string s = "letsleetcode";
    int k = 2;
    Solution solution;
    cout << solution.longestSubsequenceRepeatedK(s, k) << endl;
}
