#include "build_id.h"

#include <cstddef>
#include <cstdint>
#include <cstring>

extern "C" const std::uint8_t __start_note_gnu_build_id[];
extern "C" const std::uint8_t __stop_note_gnu_build_id[];

namespace
{

struct note_header
{
    std::uint32_t name_size;
    std::uint32_t descriptor_size;
    std::uint32_t type;
};

constexpr std::uint32_t note_type_gnu_build_id = 3;
constexpr std::uint32_t owner_size = 4;
constexpr std::uint32_t descriptor_size = 20;

// The note is a raw byte run between two linker-provided symbols, so there is no object to bind a
// reference to, and the System V gABI layout (an Elf32_Nhdr, the padded owner string, then the descriptor)
// is walked through the addresses themselves. A null return is the modeled absence of a well-formed note.
const std::uint8_t* build_id_descriptor() noexcept
{
    const std::uint8_t* note = __start_note_gnu_build_id;
    const std::size_t span = static_cast<std::size_t>(__stop_note_gnu_build_id - note);
    if(span < sizeof(note_header) + owner_size + descriptor_size)
        return nullptr;

    note_header header{};
    std::memcpy(&header, note, sizeof(header));
    if(header.type != note_type_gnu_build_id || header.name_size != owner_size || header.descriptor_size != descriptor_size)
        return nullptr;

    if(std::memcmp(note + sizeof(header), "GNU", owner_size) != 0)
        return nullptr;

    return note + sizeof(header) + owner_size;
}

}

namespace ctrlpp
{

void render_build_id(char (&text)[build_id_text_length]) noexcept
{
    const std::uint8_t* descriptor = build_id_descriptor();
    if(descriptor == nullptr)
    {
        std::memcpy(text, "unavailable", sizeof("unavailable"));
        return;
    }

    constexpr char digits[] = "0123456789abcdef";
    for(std::size_t i = 0; i < descriptor_size; ++i)
    {
        text[2 * i] = digits[descriptor[i] >> 4];
        text[2 * i + 1] = digits[descriptor[i] & 0x0Fu];
    }
    text[2 * descriptor_size] = '\0';
}

}
