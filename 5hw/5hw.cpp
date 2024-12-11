#include <bits/stdc++.h>
using namespace std;

class Solution {
public:
    int jobScheduling(vector<int> &startTime, vector<int> &endTime, vector<int> &profit) {
        size_t n = startTime.size();

        vector<vector<int> > jobs;
        // дело: {конец, начало, прибыль}
        jobs.reserve(n + 1);
        // O(n) - память

        jobs.push_back({0, 0, 0});
        // нулевая дело
        for (size_t i = 0; i < n; ++i)
            jobs.push_back({endTime[i], startTime[i], profit[i]});
        // O(n) - время

        sort(jobs.begin(), jobs.end());
        // O(n * log n) - время

        int res = 0;
        vector<int> dp;
        dp.reserve(n + 1);
        // O(n) - память
        dp.push_back(0);
        // доход от нулевой дела

        for (size_t i = 1; i < n + 1; ++i) {
            // O(n) - время

            size_t left = 1, right = i;
            while (left < right) {
                int mid = (left + right) / 2;

                if (jobs[mid][0] <= jobs[i][1])
                    left = mid + 1;
                else
                    right = mid;
            }
            // поиск последнего дела, которое заканчивается до начала текущего дела
            // O(log (right - left)) - время

            dp.push_back(max(dp[i - 1], dp[left - 1] + jobs[i][2]));
            // выбор: брать новое дело и вместе с ним все дела, которые были сделаны до его начала
            // или оставить прошлый доход
            res = max(res, dp[i]);
        }
        // O(n * log n) - время

        return res;

        // O(n * log n) - время
        // O(n) - память
    }
};
