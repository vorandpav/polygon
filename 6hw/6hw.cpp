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

        list<pair<string, vector<int> > > words;
        // список пар: {слово, массив неиспользованных букв}

        words.push_back({"", letter});
        // первое слово - пустое

        string result = "";

        while (!words.empty()) {
            int prev_size = words.size();
            // размер списка на начало шага
            auto rend = words.rend();

            for (auto word = words.rbegin(); word != rend; ++word) {
                // в порядке лексикографического возрастания

                for (int i = 0; i < 26; ++i) {
                    if (word->second[i] > 0) {
                        if (check(s, k, word->first + (char) ('a' + i))) {
                            // если новое слово встретилось k раз в строке s

                            vector<int> new_letter = word->second;
                            --new_letter[i];
                            words.push_front({word->first + (char) ('a' + i), new_letter});
                            // добавляем новое слово в начало списка, то есть самое последнее лексикографически слово будет в начале списка
                        }
                    }
                }
            }
            // O(n) - время
            result = words.front().first;
            words.resize(words.size() - prev_size);
            // удаляем слова с прошлого шага
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
