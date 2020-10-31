#include <vector>

namespace kt84 {

template <typename T>
struct Memento {
    int undo_begin;
    int undo_size;
    std::vector<T> undo_buffer;
    std::vector<T> redo_buffer;
    
    Memento(int buffer_size_ = 10) {
        init(buffer_size_);
    }
    int buffer_size() const {
        return undo_buffer.size();
    }
    void init(int buffer_size_) {
        undo_buffer.resize(buffer_size_);
        redo_buffer.clear();
        redo_buffer.reserve(buffer_size_);
        undo_begin = undo_size = 0;
    }
    void init() {
        init(buffer_size());
    }
    void store(const T& data) {
        // update undo_size or undo_begin
        if (undo_size == buffer_size())
            undo_begin = (undo_begin + 1) % buffer_size();
        else
            ++undo_size;
    
        // store current data into undo buffer
        undo_buffer[(undo_begin + undo_size - 1) % buffer_size()] = data;
        
        // clear redo buffer
        redo_buffer.clear();
    }
    bool undo(T& data) {
        if (!undo_size)
            // undo buffer is empty
            return false;
        
        // store current data into redo buffer
        redo_buffer.push_back(data);
        
        // restore data from undo buffer
        data = undo_buffer[(undo_begin + undo_size - 1) % buffer_size()];
        --undo_size;
        
        return true;
    }
    bool redo(T& data) {
        if (redo_buffer.empty())
            // redo buffer is empty
            return false;
        
        // store current data into undo buffer
        ++undo_size;
        undo_buffer[(undo_begin + undo_size - 1) % buffer_size()] = data;
        
        // restore data from redo buffer
        data = redo_buffer.back();
        redo_buffer.pop_back();
        
        return true;
    }
};

}
