#ifndef STR_UTIL_H_
#define STR_UTIL_H_

/**
 * Given a string, return an int hash for it.
 */
static inline int
hash_string(const string& s) {
	int ret = 0;
	int a = 63689;
	int b = 378551;
	for(size_t i = 0; i < s.length(); i++) {
		ret = (ret * a) + (int)s[i];
		a *= b;
	}
	return ret;
}

#endif /* STR_UTIL_H_ */
