/**
 * $Id
 */

#ifndef eman__byteorder_h__
#define eman__byteorder_h__ 1

namespace EMAN {
    
    class ByteOrder {
    public:
	static bool is_machine_big_endian();

	template <class T>
	static bool is_data_big_endian(T* small_num_addr) {
	    if (!small_num_addr) {
		return false;
	    }
	    
	    bool data_big_endian = false;
	    int typesize = sizeof(T);
	    char* d = (char*)small_num_addr;
	
	    if (is_machine_big_endian()) {
		data_big_endian = false;
		for (int i = typesize/2; i < typesize; i++) {
		    if (d[i] != 0) {
			data_big_endian = true;
			break;
		    }
		}
	    }
	    else {
		data_big_endian = true;
		for (int j = 0; j < typesize/2; j++) {
		    if (d[j] != 0) {
			data_big_endian = false;
			break;
		    }
		}
	    }
    
	    return data_big_endian;
	}
	
	
	template <class T>
	static void become_big_endian(T* data, int n = 1) {
	    if (!is_machine_big_endian()) {
		swap_bytes<T>(data, n);
	    }
	}

	template <class T>
	static void become_little_endian(T* data, int n = 1) {
	    if (is_machine_big_endian()) {
		swap_bytes<T>(data, n);
	    }
	}
	    
	
	template <class T>
	static void swap_bytes(T* data, int n = 1) {
	    unsigned char s;
	    int p = sizeof(T);
	    char* d = (char*)data;
	
	    if (p > 1) {
		for (int i = 0; i < n; i++, d += p) {
		    for (int j = 0; j < p/2; j++) {
			s = d[j];
			d[j] = d[p-1-j];
			d[p-1-j] = s;		
		    }
		}
	    }
	}
	
    private:
	static bool is_machine_endian_checked;
	static bool machine_big_endian;
    };

}

#endif
